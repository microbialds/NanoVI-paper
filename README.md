# NanoVI-paper

Benchmark scripts, result tables and figure code behind the NanoVI manuscript. NanoVI
itself, the Nextflow pipeline, is at https://github.com/microbialds/NanoVI; its
`docs/methods.md` cites the benchmark report here by section (§).

These files lived in NanoVI's repository until its 1.1.0 release (which still holds a
copy) and moved here with their history. From NanoVI 1.1.0 on, this repository is the
place where they are kept and changed.

## Layout

| Path | Contents |
|---|---|
| `benchmarks/` | The benchmark scripts (Python) and `truth/`, the expected compositions of the mock communities |
| `docs/benchmarks/methods_2026.md` | The methods benchmark report (Phase 7a, with the Phase 7b confirmation in §13); NanoVI's `docs/methods.md` cites its sections as § |
| `docs/benchmarks/m10_claims.md` | What the benchmark supports for each claim of the preprint |
| `docs/benchmarks/clinical_2026.md` | The clinical rerun (three ONT 16S cohorts) |
| `docs/benchmarks/data/` | The result tables (TSV), the source of every number in the reports and the manuscript |
| `docs/benchmarks/figures/` | The benchmark report's figures (PNG) |
| `docs/paper/` | `make_figures.py` and `style.py`, which draw the manuscript figures from `docs/benchmarks/data/` (see its README) |
| `docs/design/unknown_component.md` | Design note for the unknown-taxon component, dropped after the benchmark (a negative result, §9) |
| `docs/reviewer_response_point9.md` | Response to point 9 of the preprint's review (software engineering) |
| `nanovi/` | NanoVI, as a git submodule pinned at the v1.1.0 tag |

Paths such as `bin/`, `main.nf` or `conf/` in the reports are NanoVI's files: the
`nanovi/` checkout.

## Setup

```bash
git clone --recurse-submodules https://github.com/microbialds/NanoVI-paper
cd NanoVI-paper
conda env create -f environment.yml          # Python 3.11, matplotlib, pandas, numpy, NanoVI's pins
conda env create -f environment-tools.yml    # Badread, VSEARCH, seqkit (simulations and database variants)
conda activate nanovi-paper
```

The scripts import NanoVI's Python modules (`likelihood`, `variational_inference_bayesian`,
`em_estimator`, `build_gtdb_db`, `taxonomy`) from `nanovi/bin/`, and the clinical
rerun runs a copy of the `nanovi/` pipeline. Two variables point them elsewhere:

- `NANOVI_REPO`: another NanoVI checkout, used instead of `nanovi/` (for example
  `git worktree add ../NanoVI_34d0411 34d0411` in a NanoVI clone).
- `NANOVI_BIN`: another `bin/` only (it takes precedence over `NANOVI_REPO` for the imports).

Tools run in Docker images, which the scripts pull:

| Tool | Image | Used by |
|---|---|---|
| fastplong 0.2.2 / 0.7.0 | `quay.io/biocontainers/fastplong:0.2.2--heae3180_0`, `:0.7.0--h43da1c4_0` | `harness.py` (the pipeline's FILTER_READS) |
| minimap2 2.24 + samtools | `quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:365b17b986c1a60c1b82c6066a9345f38317b763-0` | `harness.py`, `timing.py` (ALIGN_READS) |
| minimap2 2.31 | `quay.io/biocontainers/minimap2:2.31--h118bc1c_0` | `harness.py` (the D16 tool-version arm) |
| Emu 3.6.2 | `quay.io/biocontainers/emu:3.6.2--hdfd78af_0` | `comparators.py`, `clinical.py` |
| Kraken2 2.1.3, Bracken 3.1 | `quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0`, `quay.io/biocontainers/bracken:3.1--h9948957_0` | `comparators.py` |
| BLAST 2.17.0 | `quay.io/biocontainers/blast:2.17.0--h66d330f_0` | `comparators.py` (NanoCLUST's consensus classification) |
| seqkit 2.14.0, samtools 1.24 | `quay.io/biocontainers/seqkit:2.14.0--hb192632_0`, `quay.io/biocontainers/samtools:1.24--h9dcdb79_1` | `clinical.py` |

EPI2ME wf-16s v1.6.1 (`epi2me-labs/wf-16s`) and NanoCLUST (commit 9364ddc, in the
data folder's `downloads/nanoclust`) run through Nextflow with their own containers;
NanoVI's clinical runs need Nextflow 26.04 or later. Tool versions and sources are in
§1 of `methods_2026.md`.

## Data

The data are not in this repository. They live under `~/pipelines/NanoVI_benchmarks`
on the host that ran the benchmark (`NANOVI_BENCH`, default
`~/pipelines/NanoVI_benchmarks/phase7a`; the other folders are found beside it).
Sizes in October 2026:

| Folder | Size | Contents |
|---|---|---|
| `phase7a/` | 55 GB | `downloads/` (7.7 GB: GTDB r226 SSU files, the comparators' databases, the ONT Zymo 2025.09 basecalls, SRR10391201, SRR14326315 and SRR23636353 from ENA, NanoCLUST, Zymo and *Bacillus* reference genomes; `fetch.sh` fetches them), the databases, simulations, alignments and estimates of every arm, `results/`, the conda environments and their lock files (`env_*.lock.txt`) |
| `phase7b/` | 9.3 GB | The two GTDB r226 builds (`db/`, 1.7 GB) and the Phase 7b confirmation runs (`confirm/`, 8.2 GB) |
| `phase9/` | 53 GB | Clinical cohorts: downloads (31 GB), subsampled reads (1.3 GB), runs of the four arms (21 GB), the pipeline copy that ran them (`src/`) |

Only the GTDB r226 files and the public read sets need downloading again; the
accessions and sources are in §1 and §6 of `methods_2026.md` and in
`clinical_2026.md`. The figures need none of it: `docs/paper/make_figures.py` reads
only `docs/benchmarks/data/`.

## Rerunning

| Part | Command | Writes |
|---|---|---|
| Manuscript figures | `python docs/paper/make_figures.py --out <dir>` (`--only fig3` for one) | `<name>.pdf/.svg/.png` and `<name>_data.tsv` per figure; default `--out figures/` (git-ignored) |
| Clinical tables | `python benchmarks/clinical.py analyse` | `docs/benchmarks/data/clinical_*.tsv` |
| Clinical runs | `python benchmarks/clinical.py prepare <cohort>`, then `nanovi <cohort> <gtdb\|emudb>` and `emu <cohort> [--same-sam]` (cohorts `vaginal_cl`, `vaginal_bv`, `crc`) | under `phase9/`; `prepare` copies the `nanovi/` checkout at its HEAD into `phase9/src/` |
| Database audit and variants | `python benchmarks/db_audit.py`, `python benchmarks/build_variants.py --src <exact db> --out <dir>` | `phase7a/audit/`, the variant databases |
| Simulations | `python benchmarks/simulate.py` | `phase7a/sims/` |
| Ablation grids | `python benchmarks/run_grid.py --design ofat` / `cell` | `phase7a/results/grid_*.tsv` |
| Clade exclusion, real mocks, quality filter, intervals, comparators, runtime | `clade_exclusion.py`, `real.py`, `quality_sims.py`, `intervals.py`, `run_comparators.py`, `timing.py` | `phase7a/results/` |
| Benchmark tables and figures | `python benchmarks/analyse.py <analysis>` | `docs/benchmarks/data/`, `docs/benchmarks/figures/` |
| Phase 7b confirmation | `python benchmarks/confirm_7b.py --merged <dir> --suffixes <dir>` | `docs/benchmarks/data/confirm_7b_*.tsv` |
| Database build counts | `python benchmarks/db_builds.py --merged <dir> --suffixes <dir>` | `docs/benchmarks/data/db_builds_r226.tsv` |

Each script's docstring gives its arms and options; §1 of `methods_2026.md` gives the
order of the Phase 7a stages and the compute budget.

## Which NanoVI produced each result

The results describe the NanoVI code of the commit that produced them, not the
current release. To rerun a part as it was run, check `nanovi/` out at that commit
(or point `NANOVI_REPO` at such a checkout).

| Results | NanoVI commit | Files |
|---|---|---|
| Phase 7a methods benchmark | `34d0411`, the merge of PR #13 | `methods_2026.md` §1-12, `m10_claims.md`, `docs/design/unknown_component.md`, 41 tables in `data/` (`ablation_*`, `alignment_cost_*`, `clade_exclusion_*`, `comparators_*`, `contamination_*`, `db_audit_*`, `db_variants_*`, `fig_*`, `intervals_*`, `lod_*`, `pairs_ties`, `real_*`, `timing_*`) and `figures/` |
| Phase 7b confirmation | `208349f`, the merge of PR #14 | `methods_2026.md` §13, `confirm_7b_per_sample.tsv`, `confirm_7b_summary.tsv` |
| Clinical rerun | `6659e25` (`phase9/src/COMMIT`) | `clinical_2026.md`, `clinical_*.tsv` |
| Database build counts | `6bd05f0`, the merge of PR #16 | `db_builds_r226.tsv` |
| Manuscript figures | `6bd05f0`, the merge of PR #16 | `docs/paper/` (drawn from the tables above) |

The arms whose options Phase 7b removed from NanoVI (unknown component, per-base score,
identity and aligned-length filters, `--error-profile all`, pruning off) run only with
`bin/` at `34d0411`; `harness.py` stops before any work when `bin/` lacks an option an
arm needs. With `nanovi/` at v1.1.0, `make_figures.py` and `clinical.py analyse`
reproduce the committed figure data byte for byte and the clinical tables to the last
floating-point digit (differences below 1e-15, from Python's set order).

## Citation

Cite the NanoVI preprint (Curiqueo, Fuentes-Santander and Ugalde, bioRxiv 2026,
doi:10.64898/2026.03.07.710315) and NanoVI 1.1.0; `CITATION.cff` has both.

## License

MIT, as NanoVI (`LICENSE`).
