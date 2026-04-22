# WulfRNA

WulfRNA is a packaged CLI for a focused bulk RNA-seq workflow on paired-end FASTQ input.

Current outputs:
- gene-level expected counts (aggregated from transcript quantification)
- gene-level TPM (aggregated from transcript quantification)
- final MultiQC report

Pipeline steps: FastQC (raw) → Cutadapt → FastQC (trimmed) → transcript quantification (Salmon or kallisto) → gene-level aggregation via `combined_tx2gene.tsv` → MultiQC.

## 1) Prerequisites

WulfRNA does **not** install analysis binaries for you. They must already be installed and visible in `PATH` at runtime.

Required commands:
- `bash`
- `python3`
- `fastqc`
- `cutadapt`
- `multiqc`
- `salmon` (default backend) or `kallisto` (optional backend, when selected)

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

`WORKDIR` must contain `fastq/` with strict paired-end naming:

- `sample_id_R1.fastq.gz`
- `sample_id_R2.fastq.gz`

### Reference directory

`--reference` must point to a directory containing backend-specific resources plus a shared tx2gene map:

- shared: `combined_tx2gene.tsv`
- Salmon backend: `salmon_index/`
- kallisto backend: `kallisto_index/combined_transcripts.kidx`

## 4) Run command

```bash
wulfrna run WORKDIR --reference /path/to/reference --stranded reverse --threads 16
```

Arguments:
- `WORKDIR` (positional): working directory containing `fastq/`
- `--reference PATH` (required): prepared reference directory
- `--stranded {none|forward|reverse}` (required): library strandedness
- `--threads N` (required): total threads (`N >= 1`)
- `--quantifier {salmon|kallisto}` (optional, default `salmon`)
- `--dry-run` (optional): validate tools/reference/inputs, write metadata/status, then exit
- `--genome NAME` (optional): resolve references as `<reference>/<NAME>/...`

Backward-compatibility note: legacy invocation without explicit `run` is still accepted, but `wulfrna run ...` is the intended interface.

## 5) Minimal smoke test plan

```bash
pip install -e .
wulfrna --help
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --dry-run
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --quantifier kallisto --dry-run
```

Notes:
- For dry-run to pass, required binaries must be in `PATH`, references must be complete for the selected quantifier, and `<workdir>/fastq` must contain at least one valid R1/R2 pair.

## 6) Main outputs and status markers

Expected primary outputs on full success:
- `abundance/gene_expected_counts.tsv`
- `abundance/gene_tpm.tsv`
- `multiqc/multiqc_report.html`

Status files in `WORKDIR/status/`:
- during run: `RUNNING`
- full success: `SUCCESS`, `finished_at.txt`, `summary.txt`
- dry-run success: `DRY_RUN_OK`, `finished_at.txt`, `summary.txt`
- failure: `FAILED`, `failed_step.txt` (and optionally `failed_sample.txt`), `summary.txt`
