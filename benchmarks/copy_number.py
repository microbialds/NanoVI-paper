#!/usr/bin/env python3
"""16S copy-number benchmark (Phase 15, F58, D40).

NanoVI's cell-level table divides each species' read-level abundance by its 16S
copy number (the database's copy_numbers.tsv: GTDB r226 ssu_count, the median over
the species' Complete Genome assemblies, filled from the genus or else family
median; bin/copy_numbers.py) and renormalises. This script scores the read-level
and the cell-level tables against cell-level truth, with NanoVI's own scripts
(NANOVI_REPO at a NanoVI commit with --copy-numbers, which the report names).

Database: the Phase 7a exact build (merged GTDB r226, cross-genus rule keep, D38),
rebuilt by that commit with --gtdb-metadata into <phase15>/db/merged_cn; its other
files are checked against the Phase 13 build (`check-db`).

Panels:
  real   the 12 real mocks with a cell-level truth: ONT Zymo 2025.09 MSPlus (6 reps),
         SRR10391201 and SRR23636353 (D6300/D6305), GutMock (4 reps, D6331); the
         truth is the manuals' Table 1 'Cell Number' column (benchmarks/truth/,
         cell_abundance). The Emu gut mock has none. Scores files Phase 13 stored
         (phase13/real/runs/<sample>/lik/exact.npz: fastplong 0.7.0, minimap2 2.31
         + samtools 1.24, -N 3), VB at the defaults.
  sim    seeded simulations whose reads are drawn in proportion to cells x 16S copies
         (`simulate`): Phase 7a's 20-species community, each species from one RefSeq
         'Complete Genome' assembly whose 16S copies in ssu_all are all at least
         1,400 bp (so the copies Badread draws from are the genome's ssu_count);
         no copy shared with another genus, except for a species with no genome
         that passes (S. argenteus: every complete genome has one). Cell shapes even
         (5% each) and log (Phase 7a's log weights, as cells), depths 1k and 10k,
         R9.4.1 and R10.4.1 Badread models, seeds 1-5; per sample the species' read
         counts are a multinomial draw with p ~ cells x copies. Truth: the realised
         cell shares, reads_i / copies_i renormalised. Run as the pipeline runs
         (harness.cross_run: fastplong 0.7.0, minimap2 2.31 + samtools 1.24 -N 3,
         the exact index, VB).

Arms, per sample (write_abundance_table.py --copy-numbers, one copy_numbers.tsv each):
  read     the species table (read-level), against the cell truth (and, for
           reference, against the 16S / read truth)
  db       the cell table with the database's copy numbers (D40)
  oracle   the cell table with the truth species' own copy numbers: the simulated
           genome's count, or the copy numbers Zymo's Table 1 implies (16S Only /
           Genome Copy, scaled to E. coli = 7, rounded; the database's value for the
           three D6331 species below 0.01% whose rounded shares give no ratio); other
           species keep the database's
  masked   the cell table with the truth species' own values removed (all of a
           panel's at once) and filled again by D40's rule (genus, else family median
           of the other species' values); a species whose genus and family then have
           no value drops out of the table, as designed

Metrics: L1 and Bray-Curtis (L1 / 2) at species and genus level (evaluate.compare on
the simulations, evaluate.compare_real over the species with a known share on the
mocks). `loo` also writes the database-wide leave-one-out error of D40's fill rule
(every species with its own value filled from the others) for both GTDB builds.

Outputs: <phase15>/bench/<panel>/<sample>/ (est.json, species.tsv, cell_<arm>.tsv),
<phase15>/results/copy_number.tsv (one row per sample, arm and truth),
copy_number_species.tsv, copy_number_loo.tsv, and the copy numbers used
(copy_number_values.tsv). `analyse` writes the report's tables to
docs/benchmarks/data/copy_number_*.tsv and fig_copy_number.png (its source data
under the same name).
"""
import argparse
import gzip
import json
import os
import random
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from common import BENCH, BIN, DOCS_DATA, DOCS_FIG, PHASE13, write_tsv
from evaluate import TRUTH_DIR, compare, compare_real
import copy_numbers as cn  # NanoVI bin/
import harness
import simulate

PYTHON = os.environ.get("BENCH_PYTHON", sys.executable)
PHASE15 = Path(os.environ.get("NANOVI_PHASE15", BENCH.parent / "phase15"))
DB = PHASE15 / "db" / "merged_cn"
EXACT = BENCH / "db" / "exact"
SIMS = PHASE15 / "sims"
OUT = PHASE15 / "bench"
RESULTS = PHASE15 / "results"
GENOMES = PHASE15 / "sources" / "genomes.tsv.gz"   # sources.py: GTDB r226 metadata per genome

REAL_CELL_TRUTH = {
    "zymo_MSPlus": "zymo_d6305_16s.tsv",
    "zymo_GutMock": "zymo_d6331_16s.tsv",
    "emu_zymo": "zymo_d6305_16s.tsv",
    "preprint_zymo": "zymo_d6305_16s.tsv",
}
ARMS = ["read", "db", "oracle", "masked"]
DEPTHS, CHEMISTRIES, SHAPES, SEEDS = [1000, 10000], ["r941", "r1041"], ["even", "log"], range(1, 6)
GENOME_SEED = 15


# ─── Truth and copy numbers ──────────────────────────────────────────────────

def truth_file(sample):
    for prefix, f in REAL_CELL_TRUTH.items():
        if sample.startswith(prefix):
            return TRUTH_DIR / f
    return None


def real_truth(sample, column):
    """{GTDB species: share in column}; MSPlus's spike-ins are presence only (nan), as in load_real_truth."""
    t = pd.read_csv(truth_file(sample), sep="\t")
    out = dict(zip(t["gtdb_species"], t[column]))
    if sample.startswith("zymo_MSPlus"):
        spikes = pd.read_csv(TRUTH_DIR / "zymo_msplus_spikes.tsv", sep="\t")
        out.update({s: float("nan") for s in spikes["gtdb_species"]})
    return out


def zymo_copy_numbers(truth_path, db_values):
    """Copy numbers Zymo's Table 1 implies: 16S share / genome-copy share, scaled to E. coli = 7, rounded.

    Species below 1e-4 of the 16S column (whose rounded shares give no usable ratio) take the
    database's value.
    """
    t = pd.read_csv(truth_path, sep="\t").set_index("gtdb_species")
    ratio = t["abundance"] / t["genome_copy"]
    ratio = ratio / ratio["Escherichia coli"] * 7
    return {sp: (float(round(r)) if t.loc[sp, "abundance"] >= 1e-4 else db_values[sp]) for sp, r in ratio.items()}


def db_copy_numbers():
    rows = pd.read_csv(DB / cn.COPY_NUMBER_FILE, sep="\t", dtype=str, keep_default_na=False)
    return rows


def species_values(rows):
    """species -> (copy number, n_genomes) of the species that have their own value."""
    own = rows[rows["level"] == "species"]
    return {s: (float(v), n) for s, v, n in zip(own["species"], own["copy_number"], own["n_genomes"])}


def write_variant(path, values, source):
    """copy_numbers.tsv of the database from other species values (D40's fill rule)."""
    taxonomy = cn.read_taxonomy_rows(DB / "taxonomy.tsv")
    rows = cn.fill_copy_numbers(taxonomy, values, source)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(os.devnull, "w") as quiet:
        stdout, sys.stdout = sys.stdout, quiet
        try:
            cn.write_copy_numbers(rows, path)
        finally:
            sys.stdout = stdout
    return {r["species"]: (r["copy_number"], r["level"]) for r in rows}


def copy_number_files(panel, truth_species, oracle):
    """The oracle and masked copy_numbers.tsv of a panel (one community each), written once."""
    base = species_values(db_copy_numbers())
    out = {}
    d = PHASE15 / "copy_numbers" / panel
    oracle_values = {**base, **{s: (float(v), "") for s, v in oracle.items() if v is not None}}
    masked_values = {s: v for s, v in base.items() if s not in truth_species}
    out["oracle"] = d / "oracle.tsv"
    out["masked"] = d / "masked.tsv"
    used = {}
    for arm, values in (("oracle", oracle_values), ("masked", masked_values)):
        used[arm] = write_variant(out[arm], values, f"{arm}:{panel}")
    out["db"] = DB / cn.COPY_NUMBER_FILE
    return out, used


# ─── Simulations ─────────────────────────────────────────────────────────────

def choose_complete_genomes(seed=GENOME_SEED):
    """One RefSeq Complete Genome per community species whose 16S copies are all >= 1,400 bp.

    A genome with a copy also found in another genus is skipped unless the species has no
    other (then the exception is recorded). Seeded choice among the genomes left.
    """
    audit = pd.read_csv(BENCH / "audit" / "ssu_all" / "sequences.tsv.gz", sep="\t",
                        usecols=["genome", "species", "genus", "md5", "length"])
    meta = pd.read_csv(GENOMES, sep="\t", low_memory=False).set_index("accession")
    spread = audit.drop_duplicates(["species", "md5"]).groupby("md5")["genus"].nunique()
    bad_md5 = set(spread[spread > 1].index)
    rng = random.Random(seed)
    chosen, notes = {}, {}
    for sp in simulate.community_species():
        rows = audit[audit["species"] == sp]
        per = rows.groupby("genome").agg(n=("md5", "size"), min_length=("length", "min"),
                                         cross=("md5", lambda m: bool(m.isin(bad_md5).any())))
        per["level"] = meta.reindex(per.index)["ncbi_assembly_level"]
        per["ssu_count"] = meta.reindex(per.index)["ssu_count"]
        ok = per[(per["level"] == "Complete Genome") & (per["min_length"] >= simulate.MIN_COPY_LENGTH)
                 & (per["n"] == per["ssu_count"]) & per.index.str.startswith("RS_")]
        clean = ok[~ok["cross"]]
        if len(clean) == 0:
            notes[sp] = "every complete genome has a 16S copy also found in another genus"
            clean = ok
        if len(clean) == 0:
            raise RuntimeError(f"{sp}: no RefSeq complete genome with every 16S copy >= 1,400 bp")
        chosen[sp] = rng.choice(sorted(clean.index))
    return chosen, notes


def sim_name(depth, chem, shape, seed):
    return f"cell_{shape}_{chem}_{depth // 1000}k_s{seed:02d}"


def plan(depth, chem, shape, seed, copies):
    cells = simulate.design_abundances(shape)
    species = simulate.community_species()
    w = np.array([cells[s] * copies[s] for s in species])
    rng = np.random.default_rng([15, seed, depth, CHEMISTRIES.index(chem), SHAPES.index(shape)])
    counts = rng.multinomial(depth, w / w.sum())
    return cells, dict(zip(species, counts.tolist()))


def run_simulate(threads):
    ref = SIMS / "references"
    if not (ref / "copies.tsv").exists():
        genomes, notes = choose_complete_genomes()
        df = simulate.write_references(genomes, ref)
        write_tsv(pd.DataFrame([{"species": s, "genome": g, "note": notes.get(s, "")} for s, g in genomes.items()]),
                  ref / "genomes.tsv")
    copies_df = pd.read_csv(ref / "copies.tsv", sep="\t")
    copies = copies_df.groupby("species").size().to_dict()
    species = simulate.community_species()
    samples, jobs = [], []
    for depth in DEPTHS:
        for chem in CHEMISTRIES:
            for shape in SHAPES:
                for seed in SEEDS:
                    name = sim_name(depth, chem, shape, seed)
                    out = SIMS / name
                    if (out / "reads.fastq.gz").exists():
                        continue
                    (out / "parts").mkdir(parents=True, exist_ok=True)
                    cells, counts = plan(depth, chem, shape, seed, copies)
                    samples.append((name, cells, counts))
                    for i, sp in enumerate(species):
                        if counts[sp]:
                            jobs.append((ref / f"sp{i:02d}.fasta", counts[sp], chem, 15000 + seed * 1000 + i,
                                         out / "parts" / f"sp{i:02d}.fastq"))
    print(f"{len(samples)} samples, {len(jobs)} Badread calls", flush=True)
    with ProcessPoolExecutor(threads) as pool:
        for i, _ in enumerate(pool.map(simulate.run_badread, jobs, chunksize=1), 1):
            if i % 50 == 0:
                print(f"  {i}/{len(jobs)} Badread calls done", flush=True)
    for name, cells, counts in samples:
        out = SIMS / name
        records, truth = [], []
        for i, sp in enumerate(species):
            if not counts[sp]:
                continue
            part = list(simulate.fastq_records(out / "parts" / f"sp{i:02d}.fastq"))
            if len(part) < counts[sp]:
                raise RuntimeError(f"{name} {sp}: Badread gave {len(part)} reads, {counts[sp]} needed")
            for rec in part[:counts[sp]]:
                read_id, source = rec[0][1:].split()[:2]
                truth.append((read_id, sp, source.split(",")[0]))
                records.append(rec)
        random.Random(name).shuffle(records)
        with gzip.open(out / "reads.fastq.gz", "wt", compresslevel=6) as fh:
            for rec in records:
                fh.writelines(rec)
        write_tsv(pd.DataFrame(truth, columns=["read", "species", "copy"]), out / "read_truth.tsv.gz")
        write_tsv(pd.DataFrame({"species": species, "design_cells": [cells[s] for s in species],
                                "copies": [copies[s] for s in species], "reads": [counts[s] for s in species]}),
                  out / "truth.tsv")
        for p in (out / "parts").iterdir():
            p.unlink()
        (out / "parts").rmdir()


def sim_truth(sample):
    t = pd.read_csv(SIMS / sample / "truth.tsv", sep="\t")
    reads = t["reads"] / t["reads"].sum()
    cells = (t["reads"] / t["copies"]) / (t["reads"] / t["copies"]).sum()
    return dict(zip(t["species"], reads)), dict(zip(t["species"], cells)), dict(zip(t["species"], t["copies"]))


# ─── Runs ────────────────────────────────────────────────────────────────────

def samples():
    real = sorted(p.parent.parent.name for p in (PHASE13 / "real" / "runs").glob("*/lik/exact.npz")
                  if truth_file(p.parent.parent.name))
    sims = [sim_name(d, c, s, k) for d in DEPTHS for c in CHEMISTRIES for s in SHAPES for k in SEEDS]
    return [(s, "real") for s in real] + [(s, "sim") for s in sims]


def lik_paths(sample, panel):
    if panel == "real":
        lik = PHASE13 / "real" / "runs" / sample / "lik"
        return lik / "exact.npz", lik / "exact.counts.json"
    run = OUT / panel / sample / "run"
    _, _, npz, _ = harness.cross_run(run, SIMS / sample / "reads.fastq.gz", "exact", threads=4)
    return npz, npz.with_suffix(".counts.json")


def run_sample(sample, panel, files):
    out = OUT / panel / sample
    if (out / "cell_masked.tsv").exists():
        return out
    out.mkdir(parents=True, exist_ok=True)
    npz, counts = lik_paths(sample, panel)
    est = out / "est.json"
    if not est.exists():
        subprocess.run([PYTHON, str(BIN / "estimate_abundance.py"), "--scores", str(npz), "--output", str(est) + ".tmp",
                        "--taxonomy", str(DB / "taxonomy.tsv")],
                       check=True, stdout=subprocess.DEVNULL, stderr=open(out / "est.log", "w"))
        Path(str(est) + ".tmp").rename(est)
    for arm in ARMS[1:]:
        subprocess.run([PYTHON, str(BIN / "write_abundance_table.py"), "--abundance", str(est), "--read-counts",
                        str(counts), "--taxonomy", str(DB / "taxonomy.tsv"), "--output-prefix", str(out / "species"),
                        "--counts", "--copy-numbers", str(files[arm]), "--cell-output", str(out / f"tmp_cell_{arm}")],
                       check=True, stdout=subprocess.DEVNULL, stderr=open(out / f"write_{arm}.log", "w"))
        (out / f"tmp_cell_{arm}.tsv").rename(out / f"cell_{arm}.tsv")
    return out


def _job(args):
    return run_sample(*args)


def read_table(path):
    t = pd.read_csv(path, sep="\t", comment="#")
    t = t[~t["species"].isin(["Below detection", "Unassigned"])]
    return dict(zip(t["species"], t["abundance"].astype(float)))


def n_reads(counts_path):
    return json.loads(Path(counts_path).read_text())["assigned_count"]


def score(sample, panel):
    out = OUT / panel / sample
    tables = {"read": read_table(out / "species.tsv"),
              **{arm: read_table(out / f"cell_{arm}.tsv") for arm in ARMS[1:]}}
    rows, sp_rows = [], []
    if panel == "real":
        truths = {"cell": real_truth(sample, "cell_abundance"), "16s": real_truth(sample, "abundance")}
        n = n_reads(lik_paths(sample, panel)[1])
        for arm, est in tables.items():
            for tname, truth in truths.items():
                if tname == "16s" and arm != "read":
                    continue
                m = {}
                for level in ("species", "genus"):
                    r = compare_real(est, truth, n, level)
                    m[f"{level}_l1"], m[f"{level}_bray_curtis"] = r[f"{level}_l1_known"], r[f"{level}_bray_curtis_known"]
                rows.append({"sample": sample, "panel": panel, "arm": arm, "truth": tname, **m})
        known = {k: v for k, v in truths["cell"].items() if v == v}
    else:
        reads, cells, _ = sim_truth(sample)
        truths = {"cell": cells, "16s": reads}
        for arm, est in tables.items():
            for tname, truth in truths.items():
                if tname == "16s" and arm != "read":
                    continue
                m = {}
                for level in ("species", "genus"):
                    r = compare(est, truth, level)
                    m[f"{level}_l1"], m[f"{level}_bray_curtis"] = r[f"{level}_l1"], r[f"{level}_bray_curtis"]
                rows.append({"sample": sample, "panel": panel, "arm": arm, "truth": tname, **m})
        known = cells
    # Per truth species: each table's share renormalised over the truth species, as the L1 uses them
    for arm, est in tables.items():
        s = sum(est.get(k, 0.0) for k in known)
        for sp, v in known.items():
            sp_rows.append({"sample": sample, "panel": panel, "arm": arm, "species": sp, "truth": v,
                            "estimate": est.get(sp, 0.0) / s if s else float("nan")})
    return rows, sp_rows


def run_loo():
    """Leave-one-out error of D40's fill rule over each build's species with their own value."""
    rows = []
    for build, path in (("merged", DB / cn.COPY_NUMBER_FILE),
                        ("default", PHASE15 / "cn_check" / "suffixes_keep" / cn.COPY_NUMBER_FILE),
                        ):
        if not path.exists():
            continue
        tax_path = DB / "taxonomy.tsv" if build == "merged" else PHASE13 / "db" / "suffixes_keep" / "taxonomy.tsv"
        taxonomy = cn.read_taxonomy_rows(tax_path)
        rows_cn = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        own = rows_cn[rows_cn["level"] == "species"].set_index("tax_id")["copy_number"].astype(float).to_dict()
        keys = {rank: {t["tax_id"]: cn._lineage_key(t, rank) for t in taxonomy} for rank in cn.FILL_RANKS}
        members = {rank: {} for rank in cn.FILL_RANKS}
        for rank in cn.FILL_RANKS:
            for tid in own:
                k = keys[rank][tid]
                if k is not None:
                    members[rank].setdefault(k, []).append(tid)
        for tid, value in own.items():
            filled, level = None, "none"
            for rank in cn.FILL_RANKS:
                k = keys[rank][tid]
                others = [own[o] for o in members[rank].get(k, []) if o != tid] if k is not None else []
                if others:
                    filled, level = float(np.median(others)), rank
                    break
            rows.append({"build": build, "tax_id": tid, "known": value, "filled": filled, "level": level})
    df = pd.DataFrame(rows)
    write_tsv(df, RESULTS / "copy_number_loo_species.tsv.gz")
    summ = []
    for (build, level), g in df.groupby(["build", "level"]):
        d = (g["filled"] - g["known"]).abs()
        summ.append({"build": build, "level": level, "n": len(g), "mae": d.mean(),
                     "median_rel_err": (d / g["known"]).median(), "within_1": (d <= 1).mean(),
                     "mean_signed": (g["filled"] - g["known"]).mean()})
    return write_tsv(pd.DataFrame(summ), RESULTS / "copy_number_loo.tsv")


def group_of(sample, panel):
    if panel == "real":
        if sample.startswith("zymo_GutMock"):
            return "D6331 GutMock"
        return "D6300 MSPlus" if sample.startswith("zymo_MSPlus") else "D6300 " + sample.split("_")[-1]
    shape, chem, depth = sample.split("_")[1:4]
    return f"{shape} {chem} {depth}"


def analyse():
    t = pd.read_csv(RESULTS / "copy_number.tsv", sep="\t")
    t["group"] = [group_of(s, p) for s, p in zip(t["sample"], t["panel"])]
    t["arm_truth"] = t["arm"] + " vs " + t["truth"]
    metrics = ["species_l1", "genus_l1", "species_bray_curtis"]
    # Means per panel, and per sample group
    summary = []
    for keys, by in ((["panel"], "all"), (["panel", "group"], None)):
        for k, d in t.groupby(keys + ["arm_truth"]):
            k = k if isinstance(k, tuple) else (k,)
            summary.append({"panel": k[0], "group": by or k[1], "arm_truth": k[-1], "n": len(d),
                            **{m: d[m].mean() for m in metrics}, "species_l1_sd": d["species_l1"].std()})
    summary = pd.DataFrame(summary)
    write_tsv(summary, DOCS_DATA / "copy_number_summary.tsv", float_format="%.4f")
    # Paired, sample by sample, against the cell truth: each arm against the read-level table and db
    cell = t[t["truth"] == "cell"].set_index(["panel", "sample", "arm"])["species_l1"].unstack("arm")
    paired = []
    for a, b in (("db", "read"), ("oracle", "read"), ("masked", "read"), ("db", "oracle"), ("masked", "db")):
        for panel, d in cell.groupby(level="panel"):
            diff = d[a] - d[b]
            paired.append({"panel": panel, "arm": a, "against": b, "n": len(diff), "mean_diff": diff.mean(),
                           "lower": int((diff < -1e-12).sum()), "higher": int((diff > 1e-12).sum())})
    write_tsv(pd.DataFrame(paired), DOCS_DATA / "copy_number_paired.tsv", float_format="%.4f")
    # Per species: mean |error| of each arm against the cell truth (shares over the truth species)
    sp = pd.read_csv(RESULTS / "copy_number_species.tsv", sep="\t")
    sp["abs_err"] = (sp["estimate"] - sp["truth"]).abs()
    # The copy numbers are per community: the simulations', and each Zymo truth table's
    sp["values"] = ["sim" if p == "sim" else truth_file(s).name for s, p in zip(sp["sample"], sp["panel"])]
    per = (sp.groupby(["values", "species", "arm"])
           .agg(samples=("sample", "nunique"), truth=("truth", "mean"), estimate=("estimate", "mean"),
                abs_err=("abs_err", "mean")).reset_index())
    values = pd.read_csv(RESULTS / "copy_number_values.tsv", sep="\t").rename(columns={"panel": "values"})
    per = per.merge(values[["values", "species", "db", "oracle", "masked", "masked_level"]],
                    on=["values", "species"], how="left")
    write_tsv(per, DOCS_DATA / "copy_number_species.tsv", float_format="%.4g")
    write_tsv(values, DOCS_DATA / "copy_number_values.tsv")
    write_tsv(pd.read_csv(RESULTS / "copy_number_loo.tsv", sep="\t"), DOCS_DATA / "copy_number_loo.tsv",
              float_format="%.4f")
    # The source comparison of step 1 (copy_number_sources.py): coverage, agreement, fill levels
    src = json.loads((PHASE15 / "sources" / "summary.json").read_text())
    rows = []
    for build in ("default", "merged"):
        for source, c in src[build]["coverage"].items():
            rows.append({"build": build, "species": src[build]["species"], "source": source, **c})
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "copy_number_sources_coverage.tsv")
    rows = [{"level": "genome", "pair": k.replace("agree_complete_", ""), **v}
            for k, v in src["genomes"].items() if k.startswith("agree_")]
    for build in ("default", "merged"):
        rows += [{"level": f"species ({build})", "pair": k, **v} for k, v in src[build]["agreement"].items()]
    write_tsv(pd.DataFrame(rows), DOCS_DATA / "copy_number_sources_agreement.tsv")
    g = src["genomes"]
    write_tsv(pd.DataFrame([{"assembly_level": k, "genomes": n, "mean_ssu_count": g["ssu_count_mean_by_level"][k],
                             "share_ssu_count_ge1": g["ssu_count_ge1_by_level"][k]} for k, n in g["by_level"].items()]),
              DOCS_DATA / "copy_number_sources_assembly_level.tsv")
    weighted_src = json.loads((PHASE15 / "sources" / "weighted.json").read_text())
    write_tsv(pd.DataFrame([{"samples_set": k, **{c.split(":", 1)[1]: v for c, v in w.items() if c.startswith("gtdb_complete:")},
                             "samples": w["samples"]} for k, w in weighted_src["weighted"].items()]),
              DOCS_DATA / "copy_number_sources_weighted.tsv")
    # Figure: species L1 against the cell truth per arm, real mocks by sample group and simulations by shape
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 150, "font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    arms = ["read vs 16s", "read vs cell", "db vs cell", "oracle vs cell", "masked vs cell"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
    data = []
    for ax, panel, key in ((axes[0], "real", "group"), (axes[1], "sim", None)):
        d = t[t["panel"] == panel].copy()
        if key is None:
            d["group"] = d["group"].str.split(" ").str[0] + " cells"
        groups = sorted(d["group"].unique())
        width = 0.8 / len(arms)
        for i, arm in enumerate(arms):
            m = d[d["arm_truth"] == arm].groupby("group")["species_l1"].agg(["mean", "std"]).reindex(groups)
            x = np.arange(len(groups)) + (i - (len(arms) - 1) / 2) * width
            ax.bar(x, m["mean"], width, yerr=m["std"], label=arm, capsize=1.5, error_kw={"lw": 0.5})
            data += [{"panel": panel, "group": g, "arm_truth": arm, "mean_species_l1": m.loc[g, "mean"],
                      "sd": m.loc[g, "std"]} for g in groups]
        ax.set_xticks(np.arange(len(groups)), groups, fontsize=6)
        ax.set_ylabel("species L1")
        ax.set_title("real mocks (known shares)" if panel == "real" else "simulations (reads ~ cells x copies)", fontsize=8)
    axes[1].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    DOCS_FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(DOCS_FIG / "fig_copy_number.png", bbox_inches="tight")
    write_tsv(pd.DataFrame(data), DOCS_DATA / "fig_copy_number.tsv", float_format="%.4f")
    print(summary[summary["group"] == "all"].to_string())
    print(pd.DataFrame(paired).to_string())


def check_db():
    """merged_cn must equal the Phase 13 merged keep build apart from copy_numbers.tsv."""
    ref = PHASE13 / "db" / "merged_keep"
    names = sorted(p.name for p in ref.iterdir())
    out = {}
    for name in names:
        out[name] = (ref / name).read_bytes() == (DB / name).read_bytes()
    extra = sorted({p.name for p in DB.iterdir()} - set(names))
    print(json.dumps({"same": out, "only_in_merged_cn": extra}, indent=1))
    exact_tax = (EXACT / "taxonomy.tsv").read_bytes() == (DB / "taxonomy.tsv").read_bytes()
    print("taxonomy.tsv equals the Phase 7a exact build's:", exact_tax)
    return all(out.values()) and exact_tax


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("step", choices=["check-db", "simulate", "run", "loo", "analyse"])
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=16, help="Badread processes (simulate)")
    a = ap.parse_args()
    if a.step == "check-db":
        sys.exit(0 if check_db() else 1)
    if a.step == "simulate":
        return run_simulate(a.threads)
    if a.step == "analyse":
        return analyse()
    if a.step == "loo":
        print(pd.read_csv(run_loo(), sep="\t").to_string())
        return
    harness.cross_index("exact")
    base = species_values(db_copy_numbers())
    files, used_rows = {}, []
    real_truths = {f: pd.read_csv(TRUTH_DIR / f, sep="\t") for f in set(REAL_CELL_TRUTH.values())}
    for f, t in real_truths.items():
        oracle = zymo_copy_numbers(TRUTH_DIR / f, {s: base.get(s, (None,))[0] for s in t["gtdb_species"]})
        files[f], used = copy_number_files(f.replace(".tsv", ""), set(t["gtdb_species"]), oracle)
        for sp in t["gtdb_species"]:
            used_rows.append({"panel": f, "species": sp, "db": base.get(sp, (None,))[0], "oracle": oracle[sp],
                              "masked": used["masked"][sp][0], "masked_level": used["masked"][sp][1]})
    _, _, sim_copies = sim_truth(sim_name(DEPTHS[0], CHEMISTRIES[0], SHAPES[0], SEEDS[0]))
    files["sim"], used = copy_number_files("sim", set(sim_copies), {s: float(c) for s, c in sim_copies.items()})
    for sp, c in sim_copies.items():
        used_rows.append({"panel": "sim", "species": sp, "db": base.get(sp, (None,))[0], "oracle": float(c),
                          "masked": used["masked"][sp][0], "masked_level": used["masked"][sp][1]})
    write_tsv(pd.DataFrame(used_rows), RESULTS / "copy_number_values.tsv")
    jobs = [(s, p, files["sim"] if p == "sim" else files[truth_file(s).name]) for s, p in samples()]
    with ProcessPoolExecutor(a.workers) as pool:
        for i, out in enumerate(pool.map(_job, jobs), 1):
            print(f"  {i}/{len(jobs)} {out.name}", flush=True)
    rows, sp_rows = [], []
    for s, p in samples():
        r, sr = score(s, p)
        rows += r
        sp_rows += sr
    write_tsv(pd.DataFrame(rows), RESULTS / "copy_number.tsv")
    write_tsv(pd.DataFrame(sp_rows), RESULTS / "copy_number_species.tsv")
    run_loo()
    print("done")


if __name__ == "__main__":
    main()
