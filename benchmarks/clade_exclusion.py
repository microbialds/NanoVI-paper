#!/usr/bin/env python3
"""Clade exclusion (Phase 7a, items 2b and 8; M5, D10).

Each community species is removed in turn from the exact database (every
sequence of that species; db/excl_<sp>), the index rebuilt, and the even
community's 10k-read samples (both chemistries, seeds 1-3) aligned again with the
reference setting (-N 3). On each alignment file:

  none        the pipeline as it is (Emu score, VB, floor 1e-3)
  pid<x>      alignments below x% identity skipped (Emu 3.6's --min-pid)
  len<x>      alignments with fewer than x aligned query bases skipped (--min-align-len)
  unk<δ>      the unknown-taxon component with floor at median error + δ

For each run: where the excluded species' reads go (by the responsibilities at
the final abundances: to the unknown component, to its 16S sibling in the sample,
to other community species, to species absent from the sample, unassigned or
below detection), false positives created, and the error on the species that are
in the database. Writes <bench>/results/clade_exclusion.tsv.

--database <db> (Phase 13, F56) removes each species from that database of
common.CROSS_DBS instead, and runs only the pipeline as it is ('none') with the
pipeline's tools since Phase 7b (minimap2 2.31 + samtools 1.24 on Phase 7b's
fastplong 0.7.0 reads); --seeds limits the samples. Writes under NANOVI_PHASE13
(results/clade/<db>.tsv).
"""
import argparse
import gzip
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

import harness
from common import BENCH, PHASE13, cross_db_dir, write_tsv
from evaluate import UNKNOWN, UNKNOWN_TAX_ID, compare, load_estimate, load_truth, responsibilities, taxonomy
from likelihood import read_alignment_records, score_reads
from simulate import PAIRS, community_species

SAMPLES = [f"even_{c}_10k_s{s:02d}" for c in ("r941", "r1041") for s in (1, 2, 3)]
ARMS = {
    "none": {},
    **{f"pid{x}": {"min_identity": x} for x in (85, 90, 95, 97, 98, 99)},
    **{f"len{x}": {"min_aligned_length": x} for x in (1000, 1200, 1400)},
    **{f"unk{d:g}": {"unknown_delta": d} for d in (0.01, 0.02, 0.03, 0.05)},
}
SIBLING = {a: b for p in PAIRS for a, b in (p, p[::-1])}


def db_name(sp):
    return f"excl_sp{community_species().index(sp):02d}"


def build(sp):
    db = BENCH / "db" / db_name(sp)
    if not (db / "species_taxid.fasta").exists():
        tax = taxonomy(str(BENCH / "db" / "exact"))
        tid = tax.index[tax["species"] == sp]
        excl = BENCH / "db" / f"{db_name(sp)}.exclude.txt"
        excl.write_text("\n".join(map(str, tid)) + "\n")
        subprocess.run([sys.executable, str(Path(__file__).parent / "build_variants.py"),
                        "--src", str(BENCH / "db" / "exact"), "--out", str(db), "--exclude", str(excl)], check=True)
    return harness.build_index(db, harness.index_path(db_name(sp)), threads=4)


def paths(sample, sp, arm):
    run = BENCH / "runs" / sample / "clade" / db_name(sp)
    lik = "unk" if arm.startswith("unk") else arm
    delta = ARMS[arm].get("unknown_delta")
    npz = run / (f"{lik}{delta:g}.npz" if delta else f"{lik}.npz")
    return run / "aln.bam", npz, run / f"{arm}.json"


def run_one(job):
    sample, sp, arm, threads = job
    bam, npz, est = paths(sample, sp, arm)
    harness.align(BENCH / "runs" / sample / "filtered.fastq.gz", harness.index_path(db_name(sp)), bam, 3,
                  threads=threads)
    harness.likelihood(bam, npz, **ARMS[arm])
    harness.estimate(npz, est, unknown=arm.startswith("unk"))
    return job


@lru_cache(maxsize=None)
def filtered_reads(sample, root=None):
    """Names of the reads that passed the length filter, the reads that were aligned."""
    with gzip.open((root or BENCH) / "runs" / sample / "filtered.fastq.gz", "rt") as f:
        return frozenset(line[1:].split(maxsplit=1)[0] for i, line in enumerate(f) if i % 4 == 0)


@lru_cache(maxsize=1)
def records(sample, sp):
    """The alignments of one reduced database, parsed once for its 14 arms."""
    return read_alignment_records(str(paths(sample, sp, "none")[0]))


def read_destinations(sample, sp, arm):
    """Share of the excluded species' reads by destination, from the final abundances.

    Only reads that passed the length filter count, as in the read counts that the
    in-database removal rate is computed from (analyse.py clade). The reads are scored
    again through the Python API (the same function and options as read_likelihood.py),
    which keeps their names; each read is split over its reported candidates (and the
    unknown component) in proportion to exp(score) x abundance, as the responsibilities
    at the final abundances.
    """
    bam, _, est_path = paths(sample, sp, arm)
    tax = taxonomy(str(BENCH / "db" / db_name(sp)))
    pi = {int(k): v for k, v in json.loads(est_path.read_text())["abundance"].items()}
    result = score_reads(str(bam), **ARMS[arm], records=records(sample, sp))
    log_p_rgs = result[0]
    floors = dict(zip(log_p_rgs, result[3])) if len(result) > 3 else {}
    reads = pd.read_csv(BENCH / "sims" / sample / "read_truth.tsv.gz", sep="\t")
    ours = set(reads.loc[reads["species"] == sp, "read"]) & filtered_reads(sample)
    community = set(community_species())
    dest = dict.fromkeys(["unknown", "sibling", "other community", "absent species", "below detection"], 0.0)
    for read in ours & set(log_p_rgs):
        tids, sc = log_p_rgs[read]
        cand = [(t, s) for t, s in zip(tids, sc) if pi.get(t, 0) > 0]
        if pi.get(UNKNOWN_TAX_ID, 0) > 0:
            cand.append((UNKNOWN_TAX_ID, floors[read]))
        if not cand:
            dest["below detection"] += 1
            continue
        w = np.array([s for _, s in cand]) + np.log([pi[t] for t, _ in cand])
        w = np.exp(w - w.max())
        w /= w.sum()
        for (t, _), p in zip(cand, w):
            if t == UNKNOWN_TAX_ID:
                dest["unknown"] += p
                continue
            name = tax.loc[t, "species"]
            dest["sibling" if name == SIBLING.get(sp) else
                 "other community" if name in community else "absent species"] += p
    n = len(ours)
    out = {f"to_{k.replace(' ', '_')}": v / n for k, v in dest.items()}
    out["to_unassigned"] = len(ours - set(log_p_rgs)) / n
    out["excluded_reads"] = n
    return out


def evaluate(sample, sp, arm):
    est, below = load_estimate(paths(sample, sp, arm)[2], BENCH / "db" / db_name(sp))
    truth = load_truth(BENCH / "sims" / sample)
    t = dict(zip(truth["species"], truth["realised"]))
    t_in = {k: v for k, v in t.items() if k != sp}
    tot = sum(t_in.values())
    t_in = {k: v / tot for k, v in t_in.items()}
    known = {k: v for k, v in est.items() if k != UNKNOWN}
    m = compare(known, t_in, "species")
    return {"sample": sample, "excluded": sp, "sibling_in_sample": SIBLING.get(sp, ""), "arm": arm,
            "unknown_abundance": est.get(UNKNOWN, 0.0), "below_detection": below,
            **{k.replace("species_", "indb_"): v for k, v in m.items()},
            **read_destinations(sample, sp, arm)}


# ─── Phase 13: clade exclusion on the cross-genus databases (F56) ────────────

def cross_excl_dir(db, sp):
    return PHASE13 / "db" / "excl" / db / db_name(sp)


def cross_build(db, sp):
    out = cross_excl_dir(db, sp)
    if not (out / "species_taxid.fasta").exists():
        tax = taxonomy(str(cross_db_dir(db)))
        excl = out.with_name(f"{db_name(sp)}.exclude.txt")
        excl.parent.mkdir(parents=True, exist_ok=True)
        excl.write_text("\n".join(map(str, tax.index[tax["species"] == sp])) + "\n")
        subprocess.run([sys.executable, str(Path(__file__).parent / "build_variants.py"),
                        "--src", str(cross_db_dir(db)), "--out", str(out), "--exclude", str(excl)], check=True)
    index = PHASE13 / "indexes" / "excl" / f"{db}_{db_name(sp)}.map-ont.k21.minimap2-2.31.mmi"
    return harness.build_index(out, index, threads=4, version=harness.CROSS_VERSION)


def cross_one(job):
    """Align one sample to one reduced database, score and estimate at the defaults."""
    sample, db, sp, threads = job
    run = PHASE13 / "runs" / sample
    reads = harness.cross_run(run, BENCH / "sims" / sample / "reads.fastq.gz", "exact", confirm_name=sample)[0]
    out = run / "clade" / db / db_name(sp)
    index = PHASE13 / "indexes" / "excl" / f"{db}_{db_name(sp)}.map-ont.k21.minimap2-2.31.mmi"
    harness.align(reads.resolve(), index, out / "aln.bam", 3, threads=threads, version=harness.CROSS_VERSION)
    harness.likelihood(out / "aln.bam", out / "none.npz")
    harness.estimate(out / "none.npz", out / "none.json")
    return out


def cross_evaluate(sample, db, sp):
    """Clade-exclusion metrics of one (sample, database, excluded species), as evaluate() above."""
    out = PHASE13 / "runs" / sample / "clade" / db / db_name(sp)
    tax = taxonomy(str(cross_excl_dir(db, sp)))
    est, below = load_estimate(out / "none.json", cross_excl_dir(db, sp))
    truth = load_truth(BENCH / "sims" / sample)
    t = dict(zip(truth["species"], truth["realised"]))
    t_in = {k: v for k, v in t.items() if k != sp}
    tot = sum(t_in.values())
    m = compare(est, {k: v / tot for k, v in t_in.items()}, "species")
    reads = pd.read_csv(BENCH / "sims" / sample / "read_truth.tsv.gz", sep="\t")
    ours = set(reads.loc[reads["species"] == sp, "read"]) & filtered_reads(sample, PHASE13)
    resp = responsibilities(out / "aln.bam", out / "none.json")
    community = set(community_species())
    dest = dict.fromkeys(["sibling", "other community", "absent species", "below detection"], 0.0)
    for read in ours & set(resp):
        if not resp[read]:
            dest["below detection"] += 1
        for t_id, p in resp[read].items():
            name = tax.loc[t_id, "species"]
            dest["sibling" if name == SIBLING.get(sp) else
                 "other community" if name in community else "absent species"] += p
    n = len(ours)
    return {"sample": sample, "db": db, "excluded": sp, "sibling_in_sample": SIBLING.get(sp, ""),
            "below_detection": below, **{k.replace("species_", "indb_"): v for k, v in m.items()},
            **{f"to_{k.replace(' ', '_')}": v / n for k, v in dest.items()},
            "to_unassigned": len(ours - set(resp)) / n, "excluded_reads": n,
            "align_cpu_s": harness.align_cpu(out / "aln.bam")}


def cross_main(db, species, samples, threads):
    with ProcessPoolExecutor(max(1, threads // 4)) as pool:
        list(pool.map(cross_build, [db] * len(species), species))
    jobs = [(s, db, sp, 4) for sp in species for s in samples]
    with ProcessPoolExecutor(max(1, threads // 4)) as pool:
        list(pool.map(cross_one, jobs))
    with ProcessPoolExecutor(threads) as pool:
        rows = list(pool.map(cross_evaluate, *zip(*[(s, db, sp) for s, db, sp, _ in jobs])))
    write_tsv(pd.DataFrame(rows), PHASE13 / "results" / "clade" / f"{db}.tsv")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--species", help="comma-separated subset (default: all 20)")
    ap.add_argument("--database", help="Phase 13: a database of common.CROSS_DBS (see above)")
    ap.add_argument("--seeds", default="1,2,3", help="--database: seeds of the even 10k samples (default 1,2,3)")
    args = ap.parse_args()
    species = args.species.split(",") if args.species else community_species()
    if args.database:
        seeds = {int(x) for x in args.seeds.split(",")}
        return cross_main(args.database, species, [s for s in SAMPLES if int(s[-2:]) in seeds], args.threads)
    with ProcessPoolExecutor(max(1, args.threads // 4)) as pool:
        list(pool.map(build, species))
    jobs = [(s, sp, arm, 2) for sp in species for s in SAMPLES for arm in ARMS]
    # alignments first (2 threads each), then the cheap layers reuse them
    with ProcessPoolExecutor(max(1, args.threads // 2)) as pool:
        list(pool.map(run_one, [j for j in jobs if j[2] == "none"]))
    with ProcessPoolExecutor(args.threads) as pool:
        list(pool.map(run_one, [j for j in jobs if j[2] != "none"]))
    with ProcessPoolExecutor(args.threads) as pool:
        # one chunk per (sample, species), so that a worker parses each BAM once
        rows = list(pool.map(evaluate, *zip(*[(s, sp, arm) for s, sp, arm, _ in jobs]), chunksize=len(ARMS)))
    write_tsv(pd.DataFrame(rows), BENCH / "results" / "clade_exclusion.tsv")


if __name__ == "__main__":
    main()
