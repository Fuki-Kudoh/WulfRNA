# WulfRNA
<p align="center">
  <img src="assets/wulfrna_logo.png" alt="WulfRNA logo" width="500">
</p>

WulfRNA is a packaged CLI for a focused bulk RNA-seq workflow on paired-end (default) or single-end FASTQ input.

Current outputs:
- gene-level expected counts and TPM from Salmon or kallisto
- optional coordinate-sorted STAR BAM and BAI
- STAR splice-junction and gene-count outputs
- final MultiQC report

Pipeline steps: FastQC (raw) → Cutadapt → FastQC (trimmed) →
optional STAR alignment → transcript quantification (Salmon or kallisto) →
gene-level aggregation → MultiQC.

Public-readiness note: WulfRNA v0.2.0 is a lightweight single-server
bulk RNA-seq runner intended for local or HPC workstation use.

## 1) Prerequisites

WulfRNA does **not** install analysis binaries for you. They must already be installed and visible in `PATH` at runtime.

Required commands:
- `bash`
- `python3`
- `fastqc`
- `cutadapt`
- `multiqc`
- `salmon` (default backend) or `kallisto` (optional backend, when selected)
- `STAR` and `samtools` (optional, required only when `--aligner star` is selected)

Useful standard utilities expected on Linux nodes:
- `gzip`, `zcat`, `awk`, `sed`, `grep`, `sort`, `head`, `tail`

## 2) Install the packaged CLI

Use Python 3.10+ and a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

This installs the `wulfrna` command.

## 3) Required input layout

### Work directory

`WORKDIR` must contain `fastq/` and uses an explicit input layout:

#### `paired_end` (default)

- `sample_id_R1.fastq.gz`
- `sample_id_R2.fastq.gz`

#### `single_end` (enable with `--single-end` or `--SE`)

- `sample_id.fastq.gz` **or** `sample_id_R1.fastq.gz`
- `*_R2.fastq.gz` files are not allowed in single-end mode (run fails if present).

### Reference directory

`--reference` must point to a directory containing backend-specific resources plus a shared tx2gene map:

- shared: `combined_tx2gene.tsv`
- Salmon backend: `salmon_index/`
- kallisto backend: `kallisto_index/combined_transcripts.kidx`
- STAR aligner (when `--aligner star`): `star_index/` containing non-empty `Genome`, `SA`, `SAindex`, `genomeParameters.txt`, `chrName.txt`, `chrLength.txt`, and `chrNameLength.txt`

## 4) Run command

```bash
wulfrna run WORKDIR --reference /path/to/reference --stranded reverse --threads 16
```

Arguments:
- `WORKDIR` (positional): working directory containing `fastq/`
- `--reference PATH` (required): prepared reference directory
- `--stranded {none|forward|reverse}` (required): library strandedness
  - Salmon mapping: `none -> IU`, `forward -> ISF`, `reverse -> ISR`
  - kallisto mapping: `none -> (unstranded default)`, `forward -> --fr-stranded`, `reverse -> --rf-stranded`
- `--threads N` (required): total threads (`N >= 1`)
- `--quantifier {salmon|kallisto}` (optional, default `salmon`)
- `--aligner {none|star}` (optional, default `none`): run STAR alignment from trimmed FASTQs and write isolated alignment outputs under `align/star/<sample>/`; Salmon/kallisto quantification and abundance matrices are still produced normally
- `--single-end`, `--SE` (optional): switch layout from `paired_end` to `single_end`
- `--fragment-length FLOAT` and `--fragment-sd FLOAT` (single-end kallisto only; both required and must be `> 0`)
- `--min-mapping-rate FLOAT` (optional, default `0.90`): minimum acceptable tx2gene transcript mapping rate per sample (`0.0-1.0`)
- `--dry-run` (optional): validate tools/reference/inputs, write metadata/status, then exit
- `--genome NAME` (optional): resolve references as `<reference>/<NAME>/...`
- `--no-resume` (optional): disable phase-level resume and rerun all phases
- `--force-from {fastqc_raw,cutadapt,fastqc_trimmed,align,quant,aggregate,multiqc}` (optional): force rerun from the selected phase onward

Backward-compatibility note: legacy invocation without explicit `run` is still accepted, but `wulfrna run ...` is the intended interface.

## 5) Minimal smoke test plan

```bash
pip install -e .
wulfrna --help
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --dry-run
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --quantifier kallisto --dry-run
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --aligner star --dry-run
```

Notes:
- For dry-run to pass, required binaries must be in `PATH`, references must be complete for the selected quantifier and optional aligner, and `<workdir>/fastq` must contain valid FASTQ inputs for the selected layout (`paired_end` or `single_end`).
- Single-end examples:
  - Salmon: `wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 4 --single-end --dry-run`
  - kallisto: `wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 4 --single-end --quantifier kallisto --fragment-length 200 --fragment-sd 20 --dry-run`

## 6) Main outputs and status markers

Expected primary outputs on full success:
- `abundance/gene_expected_counts.tsv`
- `abundance/gene_tpm.tsv`
- `multiqc/multiqc_report.html`
- `logs/tx2gene_mapping_stats.tsv`

Additional outputs when `--aligner star` is selected (per sample):
- `align/star/<sample>/Aligned.sortedByCoord.out.bam`
- `align/star/<sample>/Aligned.sortedByCoord.out.bam.bai`
- `align/star/<sample>/SJ.out.tab`
- `align/star/<sample>/ReadsPerGene.out.tab`
- `align/star/<sample>/Log.final.out`

STAR outputs coexist with Salmon/kallisto outputs. STAR gene counts are not aggregated into `abundance/gene_expected_counts.tsv` or `abundance/gene_tpm.tsv`.

Matrix behavior note:
- Gene-level matrices include only genes observed in the transcript quantification input (unobserved zero-only genes are not emitted).

Status files in `WORKDIR/status/`:
- during run: `RUNNING`
- full success: `SUCCESS`, `finished_at.txt`, `summary.txt`
- dry-run success: `DRY_RUN_OK`, `finished_at.txt`, `summary.txt`
- failure: `FAILED`, `failed_step.txt` (and optionally `failed_sample.txt`), `summary.txt`
- phase checkpoint markers: `status/steps/<phase>.done` for each completed phase
- run compatibility manifest: `status/manifest.json`

## 7) Automatic resume behavior

By default, `wulfrna run ...` automatically resumes at the **phase level** (not sample-level).

Phases:
- `fastqc_raw`
- `cutadapt`
- `fastqc_trimmed`
- `align` (only when `--aligner star`)
- `quant`
- `aggregate`
- `multiqc`

A phase is skipped only when:
- `status/steps/<phase>.done` exists, and
- expected outputs for that phase exist and are non-empty, and
- `status/manifest.json` is compatible with the current run configuration.

Conservative compatibility rules:
- If sample IDs changed: resume is blocked with an error.
- If input layout changed (`paired_end` vs `single_end`): automatic resume is blocked with an error.
- If `quantifier`, `reference_dir`, or `stranded` changed: `quant` and downstream phases rerun.
- If `aligner` or STAR index content changed: `align` and downstream phases rerun when STAR alignment is enabled.
- If `combined_tx2gene.tsv` fingerprint changed: `aggregate` and `multiqc` rerun.
- If only thread count changed: resume is allowed.

To fully disable resume:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 16 --no-resume
```

To rerun only aggregation/reporting after fixing `combined_tx2gene.tsv`:
- preferred: `--force-from aggregate`, or
- remove aggregate + multiqc outputs and rerun.

Example:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 16 --force-from aggregate
```
