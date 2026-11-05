"""Shared helpers for the Phase 7a benchmark scripts.

The scripts import the pipeline's own modules from NanoVI's bin/ (likelihood,
variational_inference_bayesian, em_estimator, build_gtdb_db), so the benchmark
runs the code the pipeline runs. NanoVI is the nanovi/ submodule, or the checkout
NANOVI_REPO names. Data and results live outside the repository,
under NANOVI_BENCH (default ~/pipelines/NanoVI_benchmarks/phase7a, decision D11).
NANOVI_BIN points the scripts at another bin/, such as a checkout of commit 34d0411
for the Phase 7a arms whose options Phase 7b removed (D22, D25, D29, D31).
Phase 13 (cross-genus rules, F56) writes under NANOVI_PHASE13 (default
~/pipelines/NanoVI_benchmarks/phase13) and needs NANOVI_REPO at a NanoVI commit with
--cross-genus-rule (the report names it).
"""
import gzip
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
NANOVI = Path(os.environ.get("NANOVI_REPO", REPO / "nanovi")).resolve()
BIN = Path(os.environ.get("NANOVI_BIN", NANOVI / "bin")).resolve()
if str(BIN) not in sys.path:
    sys.path.insert(0, str(BIN))

BENCH = Path(os.environ.get("NANOVI_BENCH", Path.home() / "pipelines" / "NanoVI_benchmarks" / "phase7a"))
DOCS_DATA = REPO / "docs" / "benchmarks" / "data"
DOCS_FIG = REPO / "docs" / "benchmarks" / "figures"

# Phase 13 (F56, cross-genus rules): its databases, runs and results, beside the Phase 7a data it reuses
PHASE13 = Path(os.environ.get("NANOVI_PHASE13", BENCH.parent / "phase13"))
# Phase 7b's confirmation run of the pipeline on the exact (merged) database: minimap2 2.31 + samtools 1.24
# alignments and fastplong 0.7.0 reads of the 80 simulated samples at 1k/10k and the 13 real ones
PHASE7B_CONFIRM = BENCH.parent / "phase7b" / "confirm" / "merged" / "out"
# Phase 14 (ambiguity groups, F57): its runs and results; it rescores Phase 13's stored score files
PHASE14 = Path(os.environ.get("NANOVI_PHASE14", BENCH.parent / "phase14"))
CROSS_RULES = ["remove", "best-genus", "most-copies"]
# Phase 13 databases: the Phase 7a exact build, each cross-genus rule applied to it, and c995 of each
CROSS_DBS = ["exact", *CROSS_RULES, "c995", *(f"c995_{r}" for r in CROSS_RULES)]


def cross_db_dir(db):
    """Directory of a Phase 13 database (exact is Phase 7a's merged build; c995 is rebuilt from it there)."""
    if db == "exact":
        return BENCH / "db" / "exact"
    if db == "c995":
        return BENCH / "db" / "c995"
    return PHASE13 / "db" / (db if db.startswith("c995_") else f"merged_{db}")


GTDB = BENCH / "downloads" / "gtdb_r226"
SSU_ALL = GTDB / "ssu_all_r226.fna.gz"
SSU_REPS = GTDB / "ssu_reps_r226.fna.gz"

RANKS = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]

_RAW_SPECIES = re.compile(r"s__([^;\[]*)")


def raw_gtdb_species(header):
    """GTDB species name as in the header, suffixes kept ('Listeria monocytogenes_B')."""
    m = _RAW_SPECIES.search(header)
    return m.group(1).strip() if m else ""


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_taxonomy(db_dir):
    """tax_id -> lineage dict from a NanoVI taxonomy.tsv."""
    import pandas as pd
    tax = pd.read_csv(Path(db_dir) / "taxonomy.tsv", sep="\t", dtype=str, keep_default_na=False)
    tax["tax_id"] = tax["tax_id"].astype(int)
    return tax.set_index("tax_id")


def write_tsv(df, path, **kwargs):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=kwargs.pop("index", False), **kwargs)
    return path


_RAW_RANK = {rank: re.compile(prefix + r"([^;\[]*)") for prefix, rank in
             [("d__", "superkingdom"), ("p__", "phylum"), ("c__", "class"), ("o__", "order"),
              ("f__", "family"), ("g__", "genus"), ("s__", "species")]}


def raw_gtdb_lineage(header):
    """GTDB lineage as in the header, suffixes kept, rank -> name."""
    out = {}
    for rank, rx in _RAW_RANK.items():
        m = rx.search(header)
        out[rank] = m.group(1).strip() if m else ""
    return out
