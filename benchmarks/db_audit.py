#!/usr/bin/env python3
"""Database audit (M3, M6, W9): sequences per species, 16S shared between species,
and the effect of GTDB suffix stripping, on the database build_gtdb_db.py builds.

Reads the GTDB SSU FASTA once, applying the builder's own rules (its header parser,
suffix stripping, the 900 bp minimum and exact deduplication within a species), so
the counts are those of the database the pipeline uses.

Writes, under <bench>/audit/<label>/:
  sequences.tsv.gz      one row per input sequence kept (genome, raw and stripped species, md5)
  per_species.tsv       per stripped species: raw GTDB species merged, genomes, sequences
  identical_groups.tsv  species that share at least one identical 16S sequence (connected groups)
  suffix_merges.tsv     stripped names built from more than one GTDB species, with lineage conflicts
and summary tables into docs/benchmarks/data/. With --db it summarises built
databases instead (sequences per species in each variant). With --cross-genus
(Phase 13, F56) it summarises what each cross-genus rule did, from the builds'
own reports (cross_genus_removed.tsv, cross_genus_emptied.tsv): for each build
and rule the entries and species removed, the sequences kept in one genus and
those removed from every genus (the ties a rule could not break), and the
copies lost by every species of the benchmark's truth (cross_genus_rules.tsv,
cross_genus_truth_species.tsv).
"""
import argparse
import hashlib
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from common import (BENCH, CROSS_RULES, DOCS_DATA, PHASE13, RANKS, SSU_ALL, raw_gtdb_lineage, raw_gtdb_species,
                    write_tsv)
from build_gtdb_db import fasta_records, normalize_gtdb_name, parse_ssu_header


def scan(ssu_path, min_length):
    rows = []
    for header, seq in fasta_records(str(ssu_path)):
        genome, lineage = parse_ssu_header(header)
        if lineage is None or not lineage.get("species"):
            continue
        raw = raw_gtdb_lineage(header)
        rows.append((
            genome, raw_gtdb_species(header), lineage["species"],
            *(lineage[r] for r in RANKS[1:]),
            raw["genus"], raw["family"],
            len(seq), hashlib.md5(seq.upper().encode()).hexdigest(),
        ))
    df = pd.DataFrame(rows, columns=["genome", "gtdb_species", "species", *RANKS[1:],
                                     "gtdb_genus", "gtdb_family", "length", "md5"])
    df["kept"] = df["length"] >= min_length
    return df


def union_find_groups(pairs_by_key):
    """Connected components of species linked by a shared key (an identical sequence)."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for members in pairs_by_key:
        members = list(members)
        for m in members[1:]:
            a, b = find(members[0]), find(m)
            if a != b:
                parent[b] = a
    groups = defaultdict(set)
    for x in parent:
        groups[find(x)].add(x)
    return [g for g in groups.values() if len(g) > 1]


def database_summary(db_dirs):
    """Sequences per species in built databases (exact, c995, reps), from their FASTA headers."""
    rows, hist = [], []
    for db in db_dirs:
        db = Path(db)
        counts = Counter(line[1:].split(":")[0] for line in open(db / "species_taxid.fasta") if line.startswith(">"))
        n = pd.Series(list(counts.values()))
        n_tax = sum(1 for _ in open(db / "taxonomy.tsv")) - 1
        rows.append({"database": db.name, "sequences": int(n.sum()), "species_with_sequences": len(n),
                     "species_in_taxonomy": n_tax, "mean": round(float(n.mean()), 2), "median": float(n.median()),
                     "p95": float(np.percentile(n, 95)), "p99": float(np.percentile(n, 99)), "max": int(n.max()),
                     "species_with_1": int((n == 1).sum())})
        hist.append(n.value_counts().sort_index().rename_axis("db_sequences").reset_index(name="species")
                    .assign(database=db.name))
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "db_variants_summary.tsv")
    write_tsv(pd.concat(hist)[["database", "db_sequences", "species"]], DOCS_DATA / "db_variants_seqs_per_species_hist.tsv")
    print(pd.DataFrame(rows).to_string(index=False))


def truth_species():
    """Every species the benchmark scores against: the simulated community and the real mocks' truth tables."""
    from simulate import community_species
    truth = Path(__file__).resolve().parent / "truth"
    names = set(community_species())
    for f in sorted(truth.glob("*.tsv")):
        names |= set(pd.read_csv(f, sep="\t")["gtdb_species"])
    return names


def cross_genus_summary(builds=("merged", "suffixes")):
    """What each cross-genus rule did to each GTDB r226 build (Phase 13 builds under NANOVI_PHASE13/db)."""
    rows, truth_rows = [], []
    truth = truth_species()
    for build in builds:
        keep = PHASE13 / "db" / f"{build}_keep"
        entries = Counter(line[1:].split(":")[0] for line in open(keep / "species_taxid.fasta") if line[0] == ">")
        cross = pd.read_csv(keep / "cross_genus_identical.tsv", sep="\t", dtype=str)
        tax = pd.read_csv(keep / "taxonomy.tsv", sep="\t", dtype=str).set_index("tax_id")
        groups = pd.read_csv(keep / "identical_16s.tsv", sep="\t")
        rows.append({"build": build, "rule": "keep", "entries": sum(entries.values()),
                     "species_with_sequence": len(entries), "cross_genus_sequences": len(cross),
                     "cross_genus_entries": int(cross["species_count"].astype(int).sum()),
                     "entries_removed": 0, "species_touched": 0, "species_emptied": 0,
                     "sequences_kept_in_one_genus": 0, "sequences_removed_everywhere": 0,
                     "identical_16s_species": len(groups), "identical_16s_groups": groups["group"].nunique()})
        for rule in CROSS_RULES:
            d = PHASE13 / "db" / f"{build}_{rule}"
            removed = pd.read_csv(d / "cross_genus_removed.tsv", sep="\t", dtype=str, keep_default_na=False)
            emptied = pd.read_csv(d / "cross_genus_emptied.tsv", sep="\t", dtype=str)
            groups = pd.read_csv(d / "identical_16s.tsv", sep="\t")
            n_left = sum(1 for line in open(d / "species_taxid.fasta") if line[0] == ">")
            everywhere = set(removed.loc[removed["kept_in"] == "", "md5"])
            rows.append({"build": build, "rule": rule, "entries": n_left,
                         "species_with_sequence": len(entries) - len(emptied), "cross_genus_sequences": len(cross),
                         "cross_genus_entries": int(cross["species_count"].astype(int).sum()),
                         "entries_removed": len(removed), "species_touched": removed["tax_id"].nunique(),
                         "species_emptied": len(emptied),
                         "sequences_kept_in_one_genus": len(cross) - len(everywhere),
                         "sequences_removed_everywhere": len(everywhere),
                         "identical_16s_species": len(groups), "identical_16s_groups": groups["group"].nunique()})
            lost = removed["tax_id"].value_counts()
            for t, n in lost.items():
                sp = tax.loc[t, "species"]
                if sp in truth or identical_16s_strip(sp) in truth:
                    truth_rows.append({"build": build, "rule": rule, "tax_id": t, "species": sp,
                                       "entries": entries[t], "entries_removed": int(n),
                                       "share_removed": n / entries[t], "emptied": t in set(emptied["tax_id"])})
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "cross_genus_rules.tsv")
    write_tsv(pd.DataFrame(truth_rows).sort_values(["build", "rule", "share_removed"], ascending=[True, True, False]),
              DOCS_DATA / "cross_genus_truth_species.tsv")
    print(pd.DataFrame(rows).to_string(index=False))


def identical_16s_strip(name):
    from identical_16s import strip_gtdb_suffixes
    return strip_gtdb_suffixes(name)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--ssu", default=str(SSU_ALL))
    ap.add_argument("--label", default="ssu_all")
    ap.add_argument("--min-length", type=int, default=900)
    ap.add_argument("--db", nargs="+", help="summarise these built databases instead (db_variants_summary.tsv)")
    ap.add_argument("--cross-genus", action="store_true", help="summarise the Phase 13 cross-genus builds instead")
    args = ap.parse_args()
    if args.cross_genus:
        return cross_genus_summary()
    if args.db:
        return database_summary(args.db)
    out = BENCH / "audit" / args.label

    df = scan(args.ssu, args.min_length)
    df.to_csv(out.mkdir(parents=True, exist_ok=True) or out / "sequences.tsv.gz", sep="\t", index=False)
    kept = df[df["kept"]]
    dedup = kept.drop_duplicates(["species", "md5"])  # the database: exact duplicates removed per species

    # --- Sequences per species (stripped names, as in the database) ---
    per_species = pd.DataFrame({
        "gtdb_species_merged": kept.groupby("species")["gtdb_species"].nunique(),
        "genomes": kept.groupby("species")["genome"].nunique(),
        "sequences_ge_min_length": kept.groupby("species").size(),
        "db_sequences": dedup.groupby("species").size(),
    }).fillna(0).astype(int).sort_values("db_sequences", ascending=False)
    all_species = df["species"].unique()
    n_taxonomy, n_with_seq = len(all_species), len(per_species)
    write_tsv(per_species.reset_index(), out / "per_species.tsv")

    n = per_species["db_sequences"]
    summary = [
        ("label", args.label),
        ("input sequences", len(df)),
        (f"sequences >= {args.min_length} bp", len(kept)),
        ("database sequences (exact duplicates removed per species)", len(dedup)),
        ("GTDB species in input (suffixes kept)", df["gtdb_species"].nunique()),
        ("species in taxonomy.tsv (suffixes stripped)", n_taxonomy),
        ("species with at least one database sequence", n_with_seq),
        ("species in taxonomy.tsv with no sequence", n_taxonomy - n_with_seq),
        ("sequences per species: mean", round(float(n.mean()), 2)),
        ("sequences per species: median", float(n.median())),
        ("sequences per species: 95th percentile", float(np.percentile(n, 95))),
        ("sequences per species: 99th percentile", float(np.percentile(n, 99))),
        ("sequences per species: maximum", int(n.max())),
        ("species with 1 sequence", int((n == 1).sum())),
        ("share of sequences in the top 1% of species", round(float(n.nlargest(max(1, n_with_seq // 100)).sum() / n.sum()), 4)),
    ]

    # --- Identical 16S shared between species ---
    by_md5 = dedup.groupby("md5")["species"].agg(lambda s: sorted(set(s)))
    spread = dedup.groupby("md5")[["genus", "family", "phylum"]].nunique()
    # A sequence found in species of different genera is far more likely a contaminant
    # contig in one of the assemblies than a shared 16S: groups use within-genus sharing only
    within_genus = spread["genus"] == 1
    shared_any = by_md5[by_md5.map(len) > 1]
    shared = shared_any[within_genus.loc[shared_any.index]]
    cross = spread.loc[shared_any.index]
    species_cross = set(s for m in shared_any[~within_genus.loc[shared_any.index]] for s in m)
    species_sharing = set(s for members in shared for s in members)
    # species none of whose database sequences is unique to it
    seq_species_count = by_md5.map(len)
    dedup_share = dedup.assign(n_species=dedup["md5"].map(seq_species_count))
    all_shared = dedup_share.groupby("species")["n_species"].min() > 1
    groups = union_find_groups(shared)
    group_rows = []
    for gi, g in enumerate(sorted(groups, key=len, reverse=True), 1):
        for sp in sorted(g):
            group_rows.append((f"G{gi}", len(g), sp, bool(all_shared.get(sp, False)), int(per_species.loc[sp, "db_sequences"])))
    groups_df = pd.DataFrame(group_rows, columns=["group", "group_size", "species", "no_unique_sequence", "db_sequences"])
    write_tsv(groups_df, out / "identical_groups.tsv")
    summary += [
        ("distinct database sequences", len(by_md5)),
        ("sequences identical in more than one species", len(shared_any)),
        ("  ... of which in more than one genus (likely contamination)", int((cross["genus"] > 1).sum())),
        ("  ... of which in more than one family", int((cross["family"] > 1).sum())),
        ("  ... of which in more than one phylum", int((cross["phylum"] > 1).sum())),
        ("species carrying a sequence also found in another genus", len(species_cross)),
        ("sequences identical in more than one species of the same genus", len(shared)),
        ("species sharing an identical sequence with another species of its genus", len(species_sharing)),
        ("share of species sharing an identical sequence within the genus", round(len(species_sharing) / n_with_seq, 4)),
        ("species with no sequence of their own (every variant shared, any genus)", int(all_shared.sum())),
        ("share of species with no sequence of their own", round(float(all_shared.sum()) / n_with_seq, 4)),
        ("groups of species of one genus linked by identical sequences", len(groups)),
        ("largest such group (species)", max((len(g) for g in groups), default=0)),
    ]

    # --- Suffix stripping ---
    raw_lineages = kept.drop_duplicates("gtdb_species").set_index("gtdb_species")
    raw_ranks = ["gtdb_genus", "gtdb_family", *RANKS[3:]]
    merges = kept.groupby("species")["gtdb_species"].unique()
    merges = merges[merges.map(len) > 1]
    merge_rows = []
    for sp, raws in merges.items():
        lin = raw_lineages.loc[list(raws), raw_ranks]
        conflict = [r.replace("gtdb_", "") for r in raw_ranks if lin[r].nunique() > 1]
        merge_rows.append((sp, len(raws), ";".join(sorted(raws)), ",".join(conflict),
                           conflict[-1] if conflict else "",
                           int(per_species.loc[sp, "db_sequences"])))
    merges_df = pd.DataFrame(merge_rows, columns=["species", "gtdb_species_merged", "gtdb_species",
                                                  "ranks_in_conflict", "highest_rank_in_conflict", "db_sequences"])
    write_tsv(merges_df, out / "suffix_merges.tsv")
    # genera (and higher) merged by stripping, counted over the input taxonomy
    raw_genus = df.drop_duplicates("gtdb_species")
    raw_genus_names = raw_genus.assign(raw_genus=raw_genus["gtdb_species"].str.split(" ").str[0])
    genus_merges = raw_genus_names.groupby(raw_genus_names["raw_genus"].map(normalize_gtdb_name))["raw_genus"].nunique()
    summary += [
        ("GTDB species merged by suffix stripping (stripped names with >1 source)", len(merges_df)),
        ("GTDB species absorbed into those names", int(merges_df["gtdb_species_merged"].sum())),
        ("merged names whose sources are in different GTDB genera or above (kept: the first lineage seen)",
         int((merges_df["ranks_in_conflict"] != "").sum())),
        *[(f"  ... highest conflicting rank: {r}", int((merges_df["highest_rank_in_conflict"] == r).sum())) for r in RANKS[1:]],
        ("genus names merged by suffix stripping (stripped genus with >1 source)", int((genus_merges > 1).sum())),
        ("placeholder species names (sp + digits)", int(pd.Series(all_species).str.contains(r" sp\d+$").sum())),
    ]

    summ = pd.DataFrame(summary, columns=["metric", "value"])
    write_tsv(summ, DOCS_DATA / f"db_audit_{args.label}_summary.tsv")
    write_tsv(per_species.head(20).reset_index(), DOCS_DATA / f"db_audit_{args.label}_top20.tsv")
    # source data of the distribution figure: number of species per sequence count
    write_tsv(n.value_counts().sort_index().rename_axis("db_sequences").reset_index(name="species"),
              DOCS_DATA / f"db_audit_{args.label}_seqs_per_species_hist.tsv")
    write_tsv(groups_df.groupby("group").agg(size=("group_size", "first"), species=("species", ";".join))
              .sort_values("size", ascending=False).head(30).reset_index(),
              DOCS_DATA / f"db_audit_{args.label}_identical_groups_top30.tsv")
    write_tsv(merges_df.sort_values("gtdb_species_merged", ascending=False).head(30),
              DOCS_DATA / f"db_audit_{args.label}_suffix_merges_top30.tsv")
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
