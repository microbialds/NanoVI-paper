#!/usr/bin/env python3
"""16S copy-number sources for the GTDB r226 database (Phase 15 step 1, F58, D40).

Compares GTDB r226 genome metadata (ssu_count, ncbi_ssu_count), the 16S records per
genome in ssu_all_r226.fna.gz and rrnDB 5.10 (per genome by assembly accession, and
per NCBI species through the genomes' NCBI taxonomy: the pan-taxa median of the NCBI
species most of a GTDB species' genomes carry) on the default (suffix-kept) and merged
GTDB r226 builds of Phase 13. Reports, per source, the share of the database's species
with a value, the agreement where two sources have one, the fill levels and the
leave-one-out error of filling from the genus, family and order median, and (`weighted`)
the share of real samples' assigned abundance on species with their own value.

Inputs: <phase7a>/downloads/gtdb_r226 (bac120/ar53 metadata, ssu_all) and
<phase15>/downloads/rrndb (rrnDB-5.10.tsv, rrnDB-5.10_pantaxa_stats_NCBI.tsv; the
release of 2025-04-07). Writes <phase15>/sources/: genomes.tsv.gz (per genome),
species_<build>.tsv, masked_<build>_<source>.tsv, summary.json and weighted.json.
Order-level fill is reported for the comparison; D40 fills from genus and family only.
"""
import argparse
import collections
import glob
import gzip
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

from common import BENCH, BIN, PHASE13

B = str(BENCH.parent)
GTDB = f"{BENCH}/downloads/gtdb_r226"
RRN = f"{B}/phase15/downloads/rrndb"
RANKS = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]


def strip_suffixes(name):
    import identical_16s  # NanoVI bin/
    return identical_16s.strip_gtdb_suffixes(name)


def lineage(s):
    out = {}
    for part in s.split(";"):
        p = part.strip()
        if len(p) > 3 and p[1:3] == "__":
            out[p[0]] = p[3:]
    return out


def load_metadata():
    cols = ["accession", "gtdb_taxonomy", "ssu_count", "ncbi_ssu_count", "ncbi_assembly_level",
            "ncbi_genome_category", "checkm2_completeness", "checkm2_contamination", "contig_count",
            "ncbi_taxid", "ncbi_taxonomy", "gtdb_representative"]
    frames = [pd.read_csv(f"{GTDB}/{d}_metadata_r226.tsv.gz", sep="\t", usecols=cols, low_memory=False)
              for d in ("bac120", "ar53")]
    m = pd.concat(frames, ignore_index=True)
    lin = m["gtdb_taxonomy"].map(lineage)
    m["gtdb_species"] = lin.map(lambda d: d.get("s", ""))
    m["gtdb_genus"] = lin.map(lambda d: d.get("g", ""))
    m["gtdb_family"] = lin.map(lambda d: d.get("f", ""))
    m["ncbi_species"] = m["ncbi_taxonomy"].map(lambda s: lineage(str(s)).get("s", ""))
    m["acc"] = m["accession"].str[3:]                      # GCF_000001.1
    m["acc_nov"] = m["acc"].str.split(".").str[0].str[4:]  # 000001 (GCF/GCA share the number)
    for c in ("ssu_count", "ncbi_ssu_count"):
        m[c] = pd.to_numeric(m[c], errors="coerce")
    return m


def ssu_all_counts():
    n = collections.Counter()
    with gzip.open(f"{GTDB}/ssu_all_r226.fna.gz", "rt") as fh:
        for line in fh:
            if line[0] == ">":
                n[line[1:].split("~", 1)[0]] += 1
    return pd.Series(n, name="ssu_all_records")


def load_rrndb():
    r = pd.read_csv(f"{RRN}/rrnDB-5.10.tsv", sep="\t", dtype=str)
    r = r.rename(columns={"Data source record id": "rid", "16S gene count": "rrn_16s",
                          "NCBI tax id": "rrn_taxid", "NCBI scientific name": "rrn_name"})
    r["rrn_16s"] = pd.to_numeric(r["rrn_16s"], errors="coerce")
    r["is_asm"] = r["rid"].str.match(r"^GC[AF]_\d+\.\d+$")
    r["acc_nov"] = np.where(r["is_asm"], r["rid"].str.split(".").str[0].str[4:], None)
    p = pd.read_csv(f"{RRN}/rrnDB-5.10_pantaxa_stats_NCBI.tsv", sep="\t")
    return r, p


def agreement(a, b):
    ok = a.notna() & b.notna()
    a, b = a[ok], b[ok]
    if len(a) == 0:
        return {"n": 0}
    d = a - b
    return {"n": int(len(a)), "pearson": round(float(stats.pearsonr(a, b)[0]), 4),
            "spearman": round(float(stats.spearmanr(a, b)[0]), 4),
            "mean_abs_diff": round(float(d.abs().mean()), 4), "median_abs_diff": round(float(d.abs().median()), 4),
            "equal_rounded": round(float((a.round() == b.round()).mean()), 4),
            "within_1": round(float((d.abs() <= 1).mean()), 4),
            "mean_diff_a_minus_b": round(float(d.mean()), 4)}


def fill(tax, values, ranks=("genus", "family", "order")):
    """values: Series species -> copy number. Returns (filled value, level) per species,
    filling a species without one from the median of the species values of its genus,
    then family, then order (ranks with their ancestors so a reused name is one taxon)."""
    t = tax.set_index("species")
    out = pd.DataFrame({"value": values.reindex(t.index)})
    out["level"] = np.where(out["value"].notna(), "species", "")
    for rank in ranks:
        idx = RANKS.index(rank)
        key = t[RANKS[idx:]].astype(str).agg("|".join, axis=1)
        med = out.loc[out["level"] == "species", "value"].groupby(key[out["level"] == "species"]).median()
        need = out["value"].isna()
        v = key[need].map(med)
        out.loc[need, "value"] = v
        out.loc[need & v.notna().reindex(out.index, fill_value=False), "level"] = rank
    out.loc[out["value"].isna(), "level"] = "none"
    return out


def masked_fill_error(tax, values, ranks=("genus", "family", "order")):
    """Leave one out: each species with a value is filled from the others' medians at the
    lowest rank that has another species with a value. Returns per-species rows."""
    t = tax.set_index("species").loc[values.index]
    rows = []
    keys = {r: t[RANKS[RANKS.index(r):]].astype(str).agg("|".join, axis=1) for r in ranks}
    groups = {r: collections.defaultdict(list) for r in ranks}
    for r in ranks:
        for sp, k in keys[r].items():
            groups[r][k].append(sp)
    for sp, v in values.items():
        for r in ranks:
            others = [values[o] for o in groups[r][keys[r][sp]] if o != sp]
            if others:
                f = float(np.median(others))
                rows.append((sp, v, f, r, len(others)))
                break
        else:
            rows.append((sp, v, np.nan, "none", 0))
    return pd.DataFrame(rows, columns=["species", "known", "filled", "level", "n_other"])


WEIGHTED_RANKS = ("genus", "family", "order", "class", "phylum", "superkingdom")
WEIGHTED_SETS = {
    "real_mocks": [f"{B}/phase7b/confirm/suffixes/out/sample_outputs/{p}_*_rel-abundance.tsv"
                   for p in ("emu", "preprint", "zymo")],
    "crc": [f"{B}/phase9/runs/crc/nanovi_gtdb/out/sample_outputs/*_rel-abundance.tsv"],
    "vaginal_bv": [f"{B}/phase9/runs/vaginal_bv/nanovi_gtdb/out/sample_outputs/*_rel-abundance.tsv"],
    "vaginal_cl": [f"{B}/phase9/runs/vaginal_cl/nanovi_gtdb/out/sample_outputs/*_rel-abundance.tsv"],
}


def weighted(out_dir):
    """Fill levels down to the domain (default build), and per sample the share of the assigned
    abundance whose species has its own value or takes a genus, family or order median: the 13 real
    mocks (Phase 7b's run on the default build) and Phase 9's clinical cohorts (NanoVI on GTDB r226)."""
    tax = pd.read_csv(f"{PHASE13}/db/suffixes_keep/taxonomy.tsv", sep="\t", dtype=str).fillna("")
    sp = pd.read_csv(f"{out_dir}/species_default.tsv", sep="\t", index_col=0)
    out = {"levels": {}, "weighted": {}}
    levels = {}
    for k in ("gtdb_complete", "rrndb_accession", "rrndb_ncbi_taxonomy"):
        f = fill(tax, sp[k], ranks=WEIGHTED_RANKS)
        levels[k] = f["level"]
        out["levels"][k] = {lvl: int(n) for lvl, n in f["level"].value_counts().items()}
    for name, patterns in WEIGHTED_SETS.items():
        rows = []
        for fp in sorted(f for p in patterns for f in glob.glob(p)):
            t = pd.read_csv(fp, sep="\t")
            t = t[~t["species"].isin(["Unassigned", "Below detection"])]
            tot = t["abundance"].sum()
            r = {"sample": os.path.basename(fp)}
            for k, lv in levels.items():
                lvl = t["species"].map(lv).fillna("not in db")
                for level in ("species", "genus", "family", "order"):
                    r[f"{k}:{level}"] = t.loc[lvl == level, "abundance"].sum() / tot
                r[f"{k}:above_order"] = t.loc[~lvl.isin(["species", "genus", "family", "order"]), "abundance"].sum() / tot
            rows.append(r)
        df = pd.DataFrame(rows)
        out["weighted"][name] = {"samples": len(df), **{c: round(float(df[c].mean()), 4) for c in df.columns[1:]},
                                 **{c + ":min": round(float(df[c].min()), 4) for c in df.columns[1:]
                                    if c.endswith(":species")}}
    json.dump(out, open(f"{out_dir}/weighted.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("step", nargs="?", default="compare", choices=["compare", "weighted"])
    ap.add_argument("--out", default=f"{B}/phase15/sources")
    ap.add_argument("--db", nargs=2, action="append", metavar=("NAME", "TAXONOMY"),
                    default=None)
    a = ap.parse_args()
    if a.step == "weighted":
        return weighted(a.out)
    dbs = a.db or [("default", f"{PHASE13}/db/suffixes_keep/taxonomy.tsv"),
                   ("merged", f"{PHASE13}/db/merged_keep/taxonomy.tsv")]
    os.makedirs(a.out, exist_ok=True)
    m = load_metadata()
    m["ssu_all_records"] = m["accession"].map(ssu_all_counts()).fillna(0).astype(int)
    r, pan = load_rrndb()
    rasm = r[r["is_asm"]].drop_duplicates("acc_nov")
    m = m.merge(rasm[["acc_nov", "rid", "rrn_16s"]], on="acc_nov", how="left")
    m["complete"] = m["ncbi_assembly_level"].eq("Complete Genome")
    m.to_csv(f"{a.out}/genomes.tsv.gz", sep="\t", index=False,
             columns=["accession", "gtdb_species", "gtdb_genus", "gtdb_family", "ncbi_species", "ncbi_taxid",
                      "ncbi_assembly_level", "ncbi_genome_category", "checkm2_completeness",
                      "checkm2_contamination", "contig_count", "ssu_count", "ncbi_ssu_count",
                      "ssu_all_records", "rid", "rrn_16s"])
    summary = {}
    # Per genome: what ssu_count counts
    lv = m.groupby("ncbi_assembly_level")
    summary["genomes"] = {
        "total": int(len(m)), "by_level": {k: int(v) for k, v in m["ncbi_assembly_level"].value_counts().items()},
        "ssu_count_mean_by_level": {k: round(float(g["ssu_count"].mean()), 3) for k, g in lv},
        "ssu_count_ge1_by_level": {k: round(float((g["ssu_count"] >= 1).mean()), 3) for k, g in lv},
        "ssu_count_vs_ssu_all_records_equal": round(float((m["ssu_count"] == m["ssu_all_records"]).mean()), 4),
        "ssu_count_vs_ssu_all_records_equal_complete":
            round(float((m.loc[m.complete, "ssu_count"] == m.loc[m.complete, "ssu_all_records"]).mean()), 4),
        "rrndb_records": int(len(r)), "rrndb_assembly_records": int(r["is_asm"].sum()),
        "rrndb_matched_to_gtdb": int(m["rid"].notna().sum()),
        "rrndb_matched_levels": {k: int(v) for k, v in m.loc[m["rid"].notna(), "ncbi_assembly_level"].value_counts().items()},
        "agree_complete_ssu_count_vs_rrndb": agreement(m.loc[m.complete, "ssu_count"], m.loc[m.complete, "rrn_16s"]),
        "agree_complete_ncbi_ssu_count_vs_rrndb": agreement(m.loc[m.complete, "ncbi_ssu_count"], m.loc[m.complete, "rrn_16s"]),
        "agree_complete_ssu_count_vs_ncbi_ssu_count": agreement(m.loc[m.complete, "ssu_count"], m.loc[m.complete, "ncbi_ssu_count"]),
    }
    # Species within which complete and draft genomes both exist: drafts' deficit
    comp_med = m[m.complete].groupby("gtdb_species")["ssu_count"].median()
    draft = m[~m.complete & m["ssu_count"].notna()]
    dd = draft.assign(cm=draft["gtdb_species"].map(comp_med)).dropna(subset=["cm"])
    summary["genomes"]["draft_vs_species_complete_median"] = {
        "n": int(len(dd)), "mean_draft": round(float(dd["ssu_count"].mean()), 3),
        "mean_complete_median": round(float(dd["cm"].mean()), 3),
        "share_draft_below": round(float((dd["ssu_count"] < dd["cm"]).mean()), 4),
        "share_draft_le1": round(float((dd["ssu_count"] <= 1).mean()), 4)}

    # Species values per source
    def per_species(df, col, key):
        g = df[df[col].notna() & (df[col] > 0)].groupby(key)[col]
        return g.median(), g.size()

    for name, taxpath in dbs:
        tax = pd.read_csv(taxpath, sep="\t", dtype=str).fillna("")
        merged = name == "merged"
        mm = m.copy()
        mm["sp"] = mm["gtdb_species"].map(strip_suffixes) if merged else mm["gtdb_species"]
        spset = set(tax["species"])
        src = {}
        src["gtdb_complete"], n_gc = per_species(mm[mm.complete], "ssu_count", "sp")
        src["gtdb_all"], _ = per_species(mm, "ssu_count", "sp")
        src["ssu_all_complete"], _ = per_species(mm[mm.complete], "ssu_all_records", "sp")
        src["ncbi_ssu_complete"], _ = per_species(mm[mm.complete], "ncbi_ssu_count", "sp")
        src["rrndb_accession"], n_ra = per_species(mm, "rrn_16s", "sp")
        # rrnDB through NCBI taxonomy: the NCBI species most of a GTDB species' genomes carry
        ncbi_major = (mm[mm["ncbi_species"] != ""].groupby("sp")["ncbi_species"]
                      .agg(lambda s: s.value_counts().index[0]))
        pan_sp = pan[pan["rank"] == "species"].drop_duplicates("name").set_index("name")["median"]
        src["rrndb_ncbi_taxonomy"] = ncbi_major.map(pan_sp).dropna()
        sp_tab = pd.DataFrame({k: v.reindex(sorted(spset)) for k, v in src.items()})
        sp_tab["n_complete_genomes"] = n_gc.reindex(sp_tab.index).fillna(0).astype(int)
        sp_tab["n_rrndb_genomes"] = n_ra.reindex(sp_tab.index).fillna(0).astype(int)
        sp_tab["ncbi_species_major"] = ncbi_major.reindex(sp_tab.index)
        sp_tab.index.name = "species"
        res = {"species": len(spset), "coverage": {k: {"n": int(sp_tab[k].notna().sum()),
                                                        "share": round(float(sp_tab[k].notna().mean()), 4)}
                                                    for k in src}}
        res["agreement"] = {
            "gtdb_complete_vs_rrndb_accession": agreement(sp_tab["gtdb_complete"], sp_tab["rrndb_accession"]),
            "gtdb_complete_vs_rrndb_ncbi_taxonomy": agreement(sp_tab["gtdb_complete"], sp_tab["rrndb_ncbi_taxonomy"]),
            "rrndb_accession_vs_rrndb_ncbi_taxonomy": agreement(sp_tab["rrndb_accession"], sp_tab["rrndb_ncbi_taxonomy"]),
            "gtdb_all_vs_gtdb_complete": agreement(sp_tab["gtdb_all"], sp_tab["gtdb_complete"]),
            "ssu_all_complete_vs_gtdb_complete": agreement(sp_tab["ssu_all_complete"], sp_tab["gtdb_complete"]),
        }
        # Fill rules and their masked error, per source
        res["fill"] = {}
        for k in ("gtdb_complete", "rrndb_accession", "rrndb_ncbi_taxonomy"):
            f = fill(tax, sp_tab[k])
            res["fill"][k] = {lvl: int(n) for lvl, n in f["level"].value_counts().items()}
            sp_tab[f"{k}_filled"] = f["value"]
            sp_tab[f"{k}_level"] = f["level"]
            e = masked_fill_error(tax, sp_tab[k].dropna())
            e.to_csv(f"{a.out}/masked_{name}_{k}.tsv", sep="\t", index=False)
            ee = e.dropna(subset=["filled"])
            res["fill"][k]["masked"] = {
                lvl: {"n": int(len(g)), "mae": round(float((g.filled - g.known).abs().mean()), 3),
                      "median_rel_err": round(float(((g.filled - g.known).abs() / g.known).median()), 3),
                      "within_1": round(float(((g.filled - g.known).abs() <= 1).mean()), 3)}
                for lvl, g in ee.groupby("level")}
            # Error of the rule as applied to species that would need it: species without a value
            # whose genus has others, i.e. the masked error at genus, then family
        sp_tab.to_csv(f"{a.out}/species_{name}.tsv", sep="\t", float_format="%.4g")
        summary[name] = res
    json.dump(summary, open(f"{a.out}/summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
