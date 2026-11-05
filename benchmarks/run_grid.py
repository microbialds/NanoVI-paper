#!/usr/bin/env python3
"""Ablation on identical alignments (Phase 7a, item 3: M1, M3, M4, M7, F18 -N).

Factors: database {exact, c995, reps} x -N {3, 10, 50} x error profile {primary,
all} x read score {emu, per-base} x estimator {vb, em} x output floor
(--min-abundance) {0, 1e-4, 1e-3}, plus a diagnostic with the pruning loop off
(floor 0). The reference setting is the pipeline's default: exact, -N 3, primary,
emu, vb, 1e-3.

  --design cell   every combination, on the 10k-read samples of seed 1 (both
                  chemistries and shapes): the only arms that need new alignments
  --design ofat   the reference and each of profile, score, estimator, floor and
                  pruning changed alone, on the reference alignments of every sample

Alignment against the exact database costs about 0.83 CPU-s per read (c995 0.31,
reps 0.16; -N matters little), so database and -N are varied only in the cell.
Each worker takes one sample through alignment, scores, estimates and metrics,
and writes <bench>/results/grid/<design>/<sample>.tsv (and _species.tsv.gz) as
soon as it is done; --collect merges them into grid_<design>.tsv.

  --design cross  Phase 13 (F56): the reference setting on each database of
                  common.CROSS_DBS (exact; the cross-genus rules remove, best-genus
                  and most-copies applied to it; c995 of each), with the pipeline's
                  tools since Phase 7b (fastplong 0.7.0, minimap2 2.31 + samtools
                  1.24); exact's 1k and 10k alignments are Phase 7b's pipeline run.
                  Default samples: Phase 7a's ofat set. Writes under NANOVI_PHASE13
                  (results/grid/cross/, results/grid_cross.tsv).
  --design reproduce  Phase 13's check of the Phase 7a reference: the stored exact
                  -N 3 alignments of the ofat set scored and estimated again with
                  NANOVI_REPO's bin/ (results/grid_reproduce.tsv).
"""
import argparse
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

import harness
from common import BENCH, CROSS_DBS, PHASE13, cross_db_dir, write_tsv
from evaluate import (all_metrics, cross_fp, emptied_destinations, emptied_tax_ids, load_estimate, load_truth,
                      per_species)

DBS = ["exact", "c995", "reps"]
NS = [3, 10, 50]
PROFILES = ["primary", "all"]
SCORES = ["emu", "per-base"]
ESTIMATORS = ["vb", "em"]
FLOORS = [0.0, 1e-4, 1e-3]
REFERENCE = dict(db="exact", N=3, profile="primary", score="emu", estimator="vb", floor=1e-3, prune=True)


def arms(design):
    if design == "cell":
        out = [dict(db=d, N=n, profile=p, score=s, estimator=e, floor=f, prune=True)
               for d, n, p, s, e, f in itertools.product(DBS, NS, PROFILES, SCORES, ESTIMATORS, FLOORS)]
        out += [dict(db=d, N=n, profile=p, score=s, estimator=e, floor=0.0, prune=False)
                for d, n, p, s, e in itertools.product(DBS, NS, PROFILES, SCORES, ESTIMATORS)]
        return out
    out = [dict(REFERENCE)]
    for key, values in (("profile", PROFILES), ("score", SCORES), ("estimator", ESTIMATORS), ("floor", FLOORS)):
        out += [{**REFERENCE, key: v} for v in values if v != REFERENCE[key]]
    out += [{**REFERENCE, "floor": 0.0, "prune": False}, {**REFERENCE, "estimator": "em", "floor": 0.0, "prune": False}]
    return out


def arm_label(a):
    return (f"{a['db']}_N{a['N']}_{a['profile']}_{a['score']}_{a['estimator']}_f{a['floor']:g}"
            + ("" if a["prune"] else "_noprune"))


def paths(sample, a):
    run = BENCH / "runs" / sample
    aln = run / "aln" / f"{a['db']}_N{a['N']}.bam"
    npz = run / "lik" / f"{a['db']}_N{a['N']}_{a['profile']}_{a['score']}.npz"
    est = run / "est" / f"{arm_label(a)}.json"
    return run / "filtered.fastq.gz", aln, npz, est


def sample_fields(sample):
    shape, chem, depth, seed = sample.split("_")
    return {"sample": sample, "shape": shape, "chemistry": chem, "depth": int(depth[:-1]) * 1000, "seed": int(seed[1:])}


def run_sample(sample, design, align_threads):
    """Every arm of one sample, then its metrics. Returns the per-sample result file."""
    out = BENCH / "results" / "grid" / design / f"{sample}.tsv"
    if out.exists():
        return out
    harness.filter_reads(BENCH / "sims" / sample / "reads.fastq.gz", BENCH / "runs" / sample / "filtered.fastq.gz")
    sample_arms = arms(design)
    for a in sample_arms:
        filt, aln, npz, est = paths(sample, a)
        harness.align(filt, harness.index_path(a["db"]), aln, a["N"], threads=align_threads)
        harness.likelihood(aln, npz, error_profile=a["profile"], score=a["score"])
        harness.estimate(npz, est, estimator=a["estimator"], min_abundance=a["floor"], no_pruning=not a["prune"])
    truth = load_truth(BENCH / "sims" / sample)
    rows, sp_rows = [], []
    for a in sample_arms:
        filt, aln, npz, est_path = paths(sample, a)
        est, below = load_estimate(est_path, BENCH / "db" / a["db"])
        base = {**sample_fields(sample), "arm": arm_label(a), **a}
        # Wall times under the shared machine's load: kept for reference, not performance numbers
        times = {f"wall_{k}": float(Path(str(p) + ".time").read_text()) for k, p in
                 (("align", aln), ("likelihood", npz), ("estimate", est_path))}
        rows.append({**base, **all_metrics(est, truth), "below_detection": below, **harness.read_counts(npz),
                     "bam_bytes": aln.stat().st_size, **times})
        sp_rows += [{**base, **r} for r in per_species(est, truth)]
    write_tsv(pd.DataFrame(sp_rows), out.with_name(f"{sample}_species.tsv.gz"))
    return write_tsv(pd.DataFrame(rows), out)


def cross_base(db):
    """The database a Phase 13 rule was applied to: exact, or c995 for c995 of a cleaned build."""
    return "c995" if db.startswith("c995_") else "exact"


def run_cross_sample(sample, dbs, align_threads):
    """Phase 13: every database of one simulated sample, then its metrics."""
    out = PHASE13 / "results" / "grid" / "cross" / f"{sample}.tsv"
    if out.exists():
        return out
    run = PHASE13 / "runs" / sample
    paths = {db: harness.cross_run(run, BENCH / "sims" / sample / "reads.fastq.gz", db, confirm_name=sample,
                                   threads=align_threads) for db in dbs}
    truth = load_truth(BENCH / "sims" / sample)
    expected = set(truth["species"])
    confirm_cpu = harness.confirm_align_cpu()
    rows, sp_rows = [], []
    for db in dbs:
        filt, aln, npz, est_path = paths[db]
        est, below = load_estimate(est_path, cross_db_dir(db))
        base = {**sample_fields(sample), "db": db}
        cpu = harness.align_cpu(aln) if not aln.is_symlink() else confirm_cpu.get(sample)
        moved = {}
        if db != cross_base(db) and cross_base(db) in paths:
            _, b_aln, _, b_est = paths[cross_base(db)]
            moved = emptied_destinations(b_aln, b_est, aln, est_path, emptied_tax_ids(db), expected, cross_db_dir(db))
        rows.append({**base, **all_metrics(est, truth), **cross_fp(est, expected), "below_detection": below,
                     **harness.read_counts(npz), "bam_bytes": aln.stat().st_size, "align_cpu_s": cpu,
                     "index_bytes": harness.cross_index(db).stat().st_size, **moved})
        sp_rows += [{**base, **r} for r in per_species(est, truth)]
    write_tsv(pd.DataFrame(sp_rows), out.with_name(f"{sample}_species.tsv.gz"))
    return write_tsv(pd.DataFrame(rows), out)


def run_reproduce_sample(sample):
    """Phase 13's check of the Phase 7a reference: its stored exact alignment (minimap2 2.24), scored and
    estimated again with the bin/ of NANOVI_REPO; one row with the metrics of both."""
    out = PHASE13 / "results" / "grid" / "reproduce" / f"{sample}.tsv"
    if out.exists():
        return out
    _, aln, _, ref_est = paths(sample, REFERENCE)
    npz = PHASE13 / "runs" / sample / "repro" / "exact_N3.npz"
    est = npz.with_suffix(".json")
    harness.likelihood(aln, npz)
    harness.estimate(npz, est)
    truth = load_truth(BENCH / "sims" / sample)
    row = {**sample_fields(sample)}
    for tag, path in (("7a", ref_est), ("now", est)):
        e, below = load_estimate(path, BENCH / "db" / "exact")
        row.update({f"{k}_{tag}": v for k, v in {**all_metrics(e, truth), "below_detection": below}.items()})
    return write_tsv(pd.DataFrame([row]), out)


def collect(design):
    if design == "reproduce":
        d = PHASE13 / "results" / "grid" / "reproduce"
        write_tsv(pd.concat([pd.read_csv(p, sep="\t") for p in sorted(d.glob("*.tsv"))]),
                  PHASE13 / "results" / "grid_reproduce.tsv")
        return
    if design == "cross":
        d = PHASE13 / "results" / "grid" / "cross"
        grid = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(d.glob("*.tsv"))])
        species = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(d.glob("*_species.tsv.gz"))])
        write_tsv(grid, PHASE13 / "results" / "grid_cross.tsv")
        write_tsv(species, PHASE13 / "results" / "grid_cross_species.tsv.gz")
        print(f"{len(grid)} rows from {grid['sample'].nunique()} samples -> {PHASE13}/results/grid_cross.tsv")
        return
    d = BENCH / "results" / "grid" / design
    grid = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(d.glob("*.tsv"))])
    species = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(d.glob("*_species.tsv.gz"))])
    write_tsv(grid, BENCH / "results" / f"grid_{design}.tsv")
    write_tsv(species, BENCH / "results" / f"grid_{design}_species.tsv.gz")
    print(f"{len(grid)} rows from {grid['sample'].nunique()} samples -> {BENCH}/results/grid_{design}.tsv")


def default_samples(design):
    sims = sorted(p.name for p in (BENCH / "sims").iterdir() if (p / "reads.fastq.gz").exists())
    if design == "cross":
        design = "ofat"
    if design == "cell":
        return [s for s in sims if "_10k_" in s and s.endswith("_s01")]
    # cheapest first, so that results arrive early: 1k, 10k, then 50k for seeds 1-3
    order = {"1k": 0, "10k": 1, "50k": 2}
    keep = [s for s in sims if not ("_50k_" in s and int(s[-2:]) > 3)]
    return sorted(keep, key=lambda s: (order[s.split("_")[2]], s[-3:], s))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--design", choices=["cell", "ofat", "cross", "reproduce"], required=True)
    ap.add_argument("--dbs", default=",".join(CROSS_DBS), help="--design cross: databases (default: all)")
    ap.add_argument("--samples", help="comma-separated sample names (default: by design)")
    ap.add_argument("--workers", type=int, default=4, help="samples in parallel (default 4)")
    ap.add_argument("--align-threads", type=int, default=4, help="minimap2 threads per worker (default 4)")
    ap.add_argument("--collect", action="store_true", help="only merge the per-sample results")
    args = ap.parse_args()
    if not args.collect and args.design == "reproduce":
        samples = args.samples.split(",") if args.samples else default_samples("ofat")
        with ProcessPoolExecutor(args.workers) as pool:
            list(pool.map(run_reproduce_sample, samples))
    elif not args.collect and args.design == "cross":
        dbs = args.dbs.split(",")
        for db in dbs:  # indexes first, one at a time: two workers must not build the same one
            harness.cross_index(db)
        samples = args.samples.split(",") if args.samples else default_samples("cross")
        print(f"{len(samples)} samples, {len(dbs)} databases each", flush=True)
        with ProcessPoolExecutor(args.workers) as pool:
            futures = {pool.submit(run_cross_sample, s, dbs, args.align_threads): s for s in samples}
            for f in as_completed(futures):
                try:
                    f.result()
                    print(f"done {futures[f]}", flush=True)
                except Exception as e:
                    print(f"FAILED {futures[f]}: {e!r}", flush=True)
    elif not args.collect:
        for a in arms(args.design):
            harness.check_options(dict(error_profile=a["profile"], score=a["score"]),
                                  dict(estimator=a["estimator"], min_abundance=a["floor"], no_pruning=not a["prune"]))
        samples = args.samples.split(",") if args.samples else default_samples(args.design)
        print(f"{len(samples)} samples, {len(arms(args.design))} arms each", flush=True)
        with ProcessPoolExecutor(args.workers) as pool:
            futures = {pool.submit(run_sample, s, args.design, args.align_threads): s for s in samples}
            for f in as_completed(futures):
                try:
                    f.result()
                    print(f"done {futures[f]}", flush=True)
                except Exception as e:  # keep the other samples going; the error is in the step's .log
                    print(f"FAILED {futures[f]}: {e!r}", flush=True)
    collect(args.design)


if __name__ == "__main__":
    main()
