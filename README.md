# WulfRNA

WulfRNA is a packaged CLI for a focused bulk RNA-seq workflow on paired-end FASTQ input.

Current MVP outputs:
- gene-level integer counts (featureCounts; suitable for SARTools/DESeq2)
- gene-level expected counts (RSEM)
- gene-level TPM (RSEM)
- per-sample BAM/BAI + QC artifacts
- final MultiQC report

Pipeline steps (current MVP): FastQC (raw) → Cutadapt → FastQC (trimmed) → STAR 2-pass → samtools stats/flagstat → featureCounts → RSEM → MultiQC.

## 1) Prerequisites (especially for bitwulf)

WulfRNA does **not** install analysis binaries for you. They must already be installed and visible in `PATH` at runtime.

Required commands:
- `bash`
- `python3`
- `fastqc`
- `cutadapt`
- `STAR`
- `samtools`
- `featureCounts` (from Subread)
- `rsem-calculate-expression`
- `multiqc`

Useful standard utilities expected on Linux nodes:
- `gzip`, `zcat`, `awk`, `sed`, `grep`, `sort`, `head`, `tail`

### bitwulf module + PATH setup (important)

On bitwulf, load modules **before** running `wulfrna`, then verify binaries resolve from `PATH`.

```bash
module load fastqc
module load cutadapt
module load star
module load samtools
module load subread
module load rsem
module load multiqc

which fastqc
which cutadapt
which STAR
which samtools
which featureCounts
which rsem-calculate-expression
which multiqc
```

If any `which` command prints nothing, stop and fix your module environment first.

Version checks (recommended):

```bash
fastqc --version
cutadapt --version
STAR --version
samtools --version
featureCounts -v
rsem-calculate-expression --version
multiqc --version
```

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

Example:

```text
WORKDIR/
└── fastq/
    ├── WT_1_R1.fastq.gz
    ├── WT_1_R2.fastq.gz
    ├── WT_2_R1.fastq.gz
    └── WT_2_R2.fastq.gz
```

### Reference directory

`--reference` must point to a directory containing:

- `star_index/` (STAR index directory)
- `annotation.gtf`
- `featurecounts.gtf`
- `rsem_reference.grp` and `rsem_reference.ti` (RSEM reference prefix files)

## 4) Run command (current packaged interface)

Use the CLI exactly like this pattern:

```bash
wulfrna run WORKDIR --reference /path/to/reference --stranded reverse --threads 16
```

Arguments:
- `WORKDIR` (positional): working directory containing `fastq/`
- `--reference PATH` (required): prepared reference directory
- `--stranded {none|forward|reverse}` (required): library strandedness
- `--threads N` (required): total threads (`N >= 1`)
- `--dry-run` (optional): validate tools/reference/inputs, write metadata/status, then exit
- `--genome NAME` (optional): resolve references as `<reference>/<NAME>/...`

Backward-compatibility note: legacy invocation without explicit `run` is still accepted, but `wulfrna run ...` is the intended interface.

## 5) Minimal smoke test plan

Run these in order:

```bash
pip install -e .
wulfrna --help
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --dry-run
```

Notes:
- For the dry-run smoke test to pass, required binaries must be in `PATH`, references must be complete, and `<workdir>/fastq` must contain at least one valid R1/R2 pair.
- On bitwulf, do the module loading/`which` checks first.

## 6) Main outputs and status markers

Expected primary outputs on full success:
- `counts/gene_integer_counts.tsv`
- `abundance/gene_expected_counts.tsv`
- `abundance/gene_tpm.tsv`
- `multiqc/multiqc_report.html`

Status files in `WORKDIR/status/`:
- during run: `RUNNING`
- full success: `SUCCESS`, `finished_at.txt`, `summary.txt`
- dry-run success: `DRY_RUN_OK`, `finished_at.txt`, `summary.txt`
- failure: `FAILED`, `failed_step.txt` (and optionally `failed_sample.txt`), `summary.txt`
