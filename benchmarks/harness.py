"""Steps of the abundance pipeline, run outside Nextflow so that each layer is cached
and reused across benchmark arms (Phase 7a).

The commands are those of the pipeline's modules, with the pipeline's containers:
  FILTER_READS  fastplong 0.2.2 (modules/local/filter_reads.nf, length 500-2,000,
                quality filter off unless min_quality is given)
  INDEX_DB      minimap2 2.24 -x map-ont -k 21 -d (modules/local/index_db.nf)
  ALIGN_READS   minimap2 2.24 -ax map-ont -N <N> -p 0.9 -K 500M -k 21 | samtools 1.18 view -b
  LIKELIHOOD    bin/read_likelihood.py
  RUN_VI        bin/estimate_abundance.py
The Python steps run in the benchmark environment, which has the image's package
versions (containers/requirements.txt). Version '2.31+1.24' is the pipeline since
Phase 7b (D30): minimap2 2.31 piped into samtools 1.24 in the 'alignment' label's
image, with fastplong 0.7.0 (Phase 13 runs this way); its alignments also write
<output>.cpu, the user and system CPU seconds of minimap2 and samtools. Every output is written under a temporary
name (tmp_ prefix) and renamed when complete, so an interrupted run never leaves a partial file behind,
and a step whose output exists is skipped. Each step writes <output>.time with its
wall time in seconds.
"""
import functools
import gzip
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

from common import BENCH, BIN, PHASE13, PHASE7B_CONFIRM, cross_db_dir

PYTHON = os.environ.get("BENCH_PYTHON", sys.executable)
FASTPLONG_IMAGE = {
    "0.2.2": "quay.io/biocontainers/fastplong:0.2.2--heae3180_0",
    "0.7.0": "quay.io/biocontainers/fastplong:0.7.0--h43da1c4_0",
}
ALIGN_IMAGE = {
    # minimap2 2.24 + samtools 1.18, the pipeline's 'alignment' label
    "2.24": "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0",
    # minimap2 2.31 alone (D16 arm); samtools from the 2.24 image converts its SAM
    "2.31": "quay.io/biocontainers/minimap2:2.31--h118bc1c_0",
    # minimap2 2.31 + samtools 1.24, the pipeline's 'alignment' label since Phase 7b (D30)
    "2.31+1.24": "quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:ad4c240457048f1d3cc0bd8c8fa39bf27d222eb7-0",
}
# Versions whose image has samtools, so minimap2 pipes into it as in the pipeline
PIPED = {"2.24", "2.31+1.24"}
UID = f"{os.getuid()}:{os.getgid()}"


def docker(image, cmd, mounts, cpus=None):
    # The label lets an interrupted benchmark stop its own containers
    # (docker ps -q --filter label=nanovi-phase7a | xargs docker stop)
    args = ["docker", "run", "--rm", "-u", UID, "--network", "none", "--label", "nanovi-phase7a"]
    if cpus:
        args += ["--cpus", str(cpus)]
    for m in sorted(set(mounts)):
        args += ["-v", f"{m}:{m}"]
    return args + [image, "bash", "-euo", "pipefail", "-c", cmd]


def _run(out, argv, **kwargs):
    out = Path(out)
    if out.exists():
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name("tmp_" + out.name)  # a prefix keeps the extension, which some tools read
    start = time.time()
    with open(str(out) + ".log", "w") as log:
        subprocess.run([a.replace("{OUT}", str(tmp)) for a in argv], stdout=log, stderr=subprocess.STDOUT,
                       check=True, **kwargs)
    Path(str(out) + ".time").write_text(f"{time.time() - start:.2f}\n")
    tmp.rename(out)
    return out


def filter_reads(fastq, out, min_quality=None, version="0.2.2", threads=2):
    out = Path(out)
    quality = (f"--mean_qual {min_quality} --unqualified_percent_limit 100" if min_quality is not None
               else "--disable_quality_filtering")
    report = out.with_name(out.name.replace(".fastq.gz", ""))
    cmd = (f"fastplong -i {fastq} -o {{OUT}} --thread {threads} --html {report}.html --json {report}.json "
           f"--length_required 500 --length_limit 2000 {quality}")
    return _run(out, docker(FASTPLONG_IMAGE[version], cmd, [Path(fastq).parent, out.parent]))


def build_index(db_dir, out, threads=8, version="2.24", k=21):
    out = Path(out)
    version = "2.31" if version == "2.31+1.24" else version
    cmd = f"minimap2 -x map-ont -k {k} -t {threads} -d {{OUT}} {db_dir}/species_taxid.fasta"
    return _run(out, docker(ALIGN_IMAGE[version], cmd, [Path(db_dir), out.parent]))


def align(reads, index, out, N, threads=4, version="2.24"):
    out = Path(out)
    if out.exists():  # the 2.31 path deletes its SAM, which _run would otherwise redo
        return out
    mm2 = (f"minimap2 -ax map-ont -t {threads} -N {N} -p 0.9 -K 500M -k 21 --split-prefix {out}.temp "
           f"{index} {reads}")
    if version == "2.24":
        cmd = f"{mm2} | samtools view -b -@ 2 -o {{OUT}} -"
        return _run(out, docker(ALIGN_IMAGE[version], cmd, [Path(reads).parent, Path(index).parent, out.parent]))
    if version == "2.31+1.24":
        # bash's times, not in a pipe (a subshell has no children): the second line is the
        # children's user and system CPU
        cmd = f"{mm2} | samtools view -b -@ 2 -o {{OUT}} -; times > {out}.cpu.tmp && mv {out}.cpu.tmp {out}.cpu"
        return _run(out, docker(ALIGN_IMAGE[version], cmd, [Path(reads).parent, Path(index).parent, out.parent]))
    # minimap2 2.31 writes SAM; the pipeline's samtools converts it in a second container
    sam = out.with_suffix(".sam")
    _run(sam, docker(ALIGN_IMAGE[version], f"{mm2} > {{OUT}}",
                     [Path(reads).parent, Path(index).parent, out.parent]))
    _run(out, docker(ALIGN_IMAGE["2.24"], f"samtools view -b -@ 2 -o {{OUT}} {sam}", [out.parent]))
    sam.unlink()
    return out


def subset_bam(bam, reads, out):
    """Alignments of the reads in a FASTQ, taken from a BAM of a superset of those reads.

    minimap2 maps each read independently, so this equals aligning the subset (checked:
    identical records for 250 of 500 reads on the reps index, and for 100 reads at -t 1
    against a 10k-read -t 4 alignment on the exact index). Used for filters that only
    drop reads (--min_quality).
    """
    out = Path(out)
    if out.exists():
        return out
    names = out.with_name(out.name + ".names")
    with gzip.open(reads, "rt") as fq:
        names.write_text("".join(line[1:].split()[0] + "\n" for i, line in enumerate(fq) if i % 4 == 0))
    _run(out, docker(ALIGN_IMAGE["2.24"], f"samtools view -b -N {names} -o {{OUT}} {bam}",
                     [Path(bam).parent, out.parent]))
    names.unlink()
    return out


def _likelihood_argv(bam, out, error_profile="primary", score="emu", min_identity=None, min_aligned_length=None,
                     unknown_delta=None):
    # Only options that differ from the defaults are passed: Phase 7b removed the
    # experimental ones from bin/, so those arms need bin/ at commit 34d0411 (NANOVI_BIN),
    # while the reference arm runs with the current bin/ too
    out = Path(out)
    argv = [PYTHON, str(BIN / "read_likelihood.py"), "--alignments", str(bam), "--scores", "{OUT}",
            "--counts", str(out.with_suffix(".counts.json"))]
    for flag, value in (("--error-profile", None if error_profile == "primary" else error_profile),
                        ("--score", None if score == "emu" else score),
                        ("--min-identity", min_identity), ("--min-aligned-length", min_aligned_length),
                        ("--unknown-delta", unknown_delta)):
        if value is not None:
            argv += [flag, str(value)]
    return argv


def likelihood(bam, out, **options):
    return _run(out, _likelihood_argv(bam, out, **options))


def _estimate_argv(npz, out, estimator="vb", min_abundance=0.001, no_pruning=False, unknown=False):
    argv = [PYTHON, str(BIN / "estimate_abundance.py"), "--scores", str(npz), "--output", "{OUT}",
            "--estimator", estimator, "--min-abundance", str(min_abundance)]
    if no_pruning:
        argv.append("--no-pruning")
    if unknown:
        argv.append("--unknown-component")
    return argv


@functools.lru_cache(maxsize=None)
def _help(script):
    return subprocess.run([PYTHON, str(BIN / script), "--help"], capture_output=True, text=True).stdout


def check_options(likelihood_options, estimate_options):
    """Stop before any work when bin/ lacks an option an arm passes.

    The arms whose options Phase 7b removed (D22, D25, D29, D31) need NANOVI_BIN set
    to bin/ of commit 34d0411; without it they would fail only when their turn came.
    """
    for script, argv in (("read_likelihood.py", _likelihood_argv("x.bam", "x.npz", **likelihood_options)),
                         ("estimate_abundance.py", _estimate_argv("x.npz", "x.json", **estimate_options))):
        # '--score' must not match '--scores'
        missing = [a for a in argv
                   if a.startswith("--") and not re.search(re.escape(a) + r"(?![\w-])", _help(script))]
        if missing:
            sys.exit(f"{BIN / script} has no {', '.join(missing)}: this arm needs NANOVI_BIN set to bin/ "
                     "of commit 34d0411 (Phase 7b removed the option)")


def estimate(npz, out, **options):
    argv = _estimate_argv(npz, out, **options)
    if not Path(out).exists() and read_counts(npz)["assigned_count"] == 0:
        # a filter that keeps no reads: the pipeline stops here; record an empty estimate
        # (tmp_ and rename as in _run, with a .time file as every other step has)
        tmp = Path(out).with_name("tmp_" + Path(out).name)
        tmp.write_text(json.dumps({"abundance": {}, "below_detection_count": 0, "read_summary": {}}))
        Path(str(out) + ".time").write_text("0.00\n")
        tmp.rename(out)
        return Path(out)
    return _run(out, argv, env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                                "MKL_NUM_THREADS": "1"})


def read_counts(npz):
    return json.loads(Path(npz).with_suffix(".counts.json").read_text())


def index_path(db, version="2.24", k=21, root=None):
    # one index per minimap2 version: the pipeline's image builds the same index as minimap2 2.31's
    version = version.split("+")[0]
    suffix = "" if version == "2.24" else f".minimap2-{version}"
    return (root or BENCH) / "indexes" / f"{db}.map-ont.k{k}{suffix}.mmi"


def align_cpu(bam):
    """User + system CPU seconds of an alignment run with version '2.31+1.24' (None for others)."""
    path = Path(str(bam) + ".cpu")
    if not path.exists():
        return None
    # second line of bash's times: children's user and system CPU, as 'XmY.YYYs XmY.YYYs'
    return sum(60 * float(m) + float(sec) for m, sec in re.findall(r"(\d+)m([\d.]+)s",
                                                                   path.read_text().splitlines()[1]))


# ─── Phase 13: cross-genus database arms (F56) ──────────────────────────────

CROSS_VERSION = "2.31+1.24"


def cross_index(db):
    """minimap2 2.31 index of a Phase 13 database; exact's is the Phase 7a/7b one."""
    root = BENCH if db == "exact" else PHASE13
    path = index_path(db, CROSS_VERSION, root=root)
    if not path.exists():
        build_index(cross_db_dir(db), path, threads=8, version=CROSS_VERSION)
    return path


def _link(src, dst):
    dst = Path(dst)
    if not dst.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.symlink_to(src)
    return dst


def cross_run(run_dir, raw_reads, db, confirm_name=None, threads=4):
    """One Phase 13 arm as the pipeline runs it: fastplong 0.7.0, minimap2 2.31 + samtools 1.24 at
    -N 3, the Emu score, VB at the defaults. Returns (filtered reads, BAM, npz, estimate).

    Where Phase 7b's pipeline run of the exact database has the sample (confirm_name), its
    filtered reads are used for every database and its alignment for exact, so the exact arm
    is that run.
    """
    run = Path(run_dir)
    filt, aln, npz, est = (run / "filtered.fastq.gz", run / "aln" / f"{db}.bam",
                           run / "lik" / f"{db}.npz", run / "est" / f"{db}.json")
    confirm_reads = PHASE7B_CONFIRM / "filtered" / f"{confirm_name}_filtered.fastq.gz"
    if confirm_name and confirm_reads.exists():
        _link(confirm_reads, filt)
        if db == "exact":
            _link(PHASE7B_CONFIRM / "bam" / f"{confirm_name}_vi_alignments.bam", aln)
    else:
        filter_reads(raw_reads, filt, version="0.7.0")
    # a linked file's own directory is what the container must mount
    align(filt.resolve(), cross_index(db), aln, 3, threads=threads, version=CROSS_VERSION)
    likelihood(aln, npz)
    estimate(npz, est)
    return filt, aln, npz, est


@functools.lru_cache(maxsize=None)
def confirm_align_cpu():
    """sample -> CPU seconds of ALIGN_READS in Phase 7b's pipeline run (realtime x %cpu, its trace)."""
    import csv
    unit = {"h": 3600, "m": 60, "s": 1, "ms": 0.001}
    out = {}
    with open(PHASE7B_CONFIRM / "pipeline_info" / "trace.txt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if "ALIGN_READS" not in row["name"]:
                continue
            secs = sum(float(v) * unit[u] for v, u in re.findall(r"([\d.]+)(ms|h|m|s)", row["realtime"]))
            # the trace names tasks by number: the BAM in the task's work directory names the sample
            work = next((PHASE7B_CONFIRM.parent / "work" / row["hash"].split("/")[0]).glob(
                row["hash"].split("/")[1] + "*"))
            sample = next(work.glob("*_vi_alignments.bam")).name.replace("_vi_alignments.bam", "")
            out[sample] = secs * float(row["%cpu"].rstrip("%")) / 100
    return out
