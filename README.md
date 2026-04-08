# WulfRNA

WulfRNA is a bulk RNA-seq pipeline for **paired-end** FASTQ files on a single-server environment.

It is designed for a simple lab workflow:

1. Put paired-end FASTQ files in `fastq/`
2. Run one command
3. Get:
   - gene-level integer counts for **SARTools / DESeq2**
   - gene-level **expected counts**
   - gene-level **TPM**
   - BAM files
   - a final **MultiQC** report

Current pipeline scope:

- FastQC (raw)
- Cutadapt
- FastQC (trimmed)
- STAR 2-pass
- samtools stats
- featureCounts
- RSEM
- MultiQC

---

## 1. What must be installed before running WulfRNA

WulfRNA does **not** install analysis tools for you.  
Before running the pipeline, you must ensure that the following commands are installed and available in your `PATH`.

### Required commands

- `bash`
- `python3`
- `fastqc`
- `cutadapt`
- `STAR`
- `samtools`
- `featureCounts`
- `rsem-calculate-expression`
- `multiqc`

### Strongly recommended standard utilities

These are usually already present on a Linux server, but WulfRNA should assume they are available:

- `gzip`
- `zcat`
- `awk`
- `sed`
- `grep`
- `sort`
- `head`
- `tail`

### Notes on tool provenance

- `featureCounts` is usually provided by the **Subread** package.
- `RSEM` should be installed with its normal command-line tools.
- `STAR` requires a prebuilt genome index.
- `RSEM` requires a prebuilt RSEM reference.
- `MultiQC` runs only once at the end, after all sample-level steps complete successfully.

---

## 2. bitwulf-specific setup

On bitwulf, tools are managed with **environment modules**.  
That means WulfRNA should **not** hardcode absolute executable paths into the code.  
Instead, load the required modules first, then confirm the commands are visible in `PATH`.

### Example module loading

Adjust the module names below to match your actual bitwulf module names.

```bash
module load fastqc
module load cutadapt
module load star
module load samtools
module load subread
module load rsem
module load multiqc
````

Depending on your environment, the real module names may differ. For example, a tool may live under a versioned module such as `samtools/1.23` rather than `samtools`.

### Confirm that commands are visible

Run:

```bash
which fastqc
which cutadapt
which STAR
which samtools
which featureCounts
which rsem-calculate-expression
which multiqc
```

Each command should print a valid path.

You can also confirm versions:

```bash
fastqc --version
cutadapt --version
STAR --version
samtools --version
featureCounts -v
rsem-calculate-expression --version
multiqc --version
```

If any command is missing, the pipeline is **not ready to run**.

---

## 3. Python environment for WulfRNA itself

WulfRNA is expected to run with `python3`.

Recommended:

* Python 3.10 or newer
* a dedicated virtual environment for the pipeline code

Example:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the package in editable mode:

```bash
pip install -e .
```

This installs the `wulfrna` command.
You can also run the CLI module directly:

```bash
python -m wulfrna.cli run /path/to/workdir --reference /path/to/reference_dir --stranded reverse --threads 16
```

---

## 4. Reference files required before running

WulfRNA expects a reference directory prepared **before** execution.

The reference directory must contain at least:

* a **STAR genome index**
* a **gene annotation GTF**
* a **featureCounts-compatible annotation GTF**
* an **RSEM reference prefix**

Example conceptual layout:

```text
/reference_dir/
├── star_index/
│   ├── Genome
│   ├── SA
│   ├── SAindex
│   └── ...
├── annotation/
│   └── genes.gtf
└── rsem/
    ├── reference.grp
    ├── reference.ti
    ├── reference.transcripts.fa
    └── ...
```

The exact filenames may vary depending on implementation, but the pipeline must be able to resolve the equivalent resources.

If the reference directory is incomplete, the pipeline should fail **before** any sample is processed.

---

## 5. Input FASTQ naming rules

The working directory must contain a `fastq/` directory.

Inside `fastq/`, paired-end FASTQ files must follow this naming convention:

* `sample_id_R1.fastq.gz`
* `sample_id_R2.fastq.gz`

Examples:

```text
fastq/
├── WT_1_R1.fastq.gz
├── WT_1_R2.fastq.gz
├── WT_2_R1.fastq.gz
├── WT_2_R2.fastq.gz
├── KO_1_R1.fastq.gz
└── KO_1_R2.fastq.gz
```

### Important rules

* Only gzipped FASTQ is supported: `*.fastq.gz`
* The suffix must be exactly `_R1.fastq.gz` and `_R2.fastq.gz`
* Each `R1` file must have a matching `R2`
* Sample IDs are derived from the prefix before `_R1.fastq.gz` or `_R2.fastq.gz`

If the input naming is inconsistent, the pipeline should fail immediately.

---

## 6. Working directory layout

The user prepares the working directory like this:

```text
workdir/
└── fastq/
    ├── sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz
    ├── sample2_R1.fastq.gz
    └── sample2_R2.fastq.gz
```

During execution, WulfRNA will create additional directories such as:

```text
workdir/
├── fastq/
├── trimmed/
├── align/
├── qc/
├── counts/
├── abundance/
├── multiqc/
├── logs/
├── metadata/
├── status/
└── samples.tsv
```

---

## 7. How to run the pipeline

The intended interface is:

```bash
wulfrna run /path/to/workdir \
  --reference /path/to/reference_dir \
  --stranded reverse \
  --threads 16
```

### Arguments

* `/path/to/workdir`
  Working directory containing `fastq/`

* `--reference /path/to/reference_dir`
  Directory containing STAR index, annotation GTF, and RSEM reference

* `--stranded {none|forward|reverse}`
  Library strandedness

* `--threads N`
  Number of CPU threads to use

### Example

```bash
wulfrna run /data/fuki/projects/test_run \
  --reference /data/fuki/references/mm10_M21 \
  --stranded reverse \
  --threads 16
```

If the code has not yet been packaged into a CLI, the equivalent development form may be:

```bash
python3 -m wulfrna.cli run /data/fuki/projects/test_run \
  --reference /data/fuki/references/mm10_M21 \
  --stranded reverse \
  --threads 16
```

Use whichever entrypoint is actually implemented in the repository.

---

## 8. What the pipeline does

For all paired-end samples detected in `fastq/`, WulfRNA will:

1. run **FastQC** on raw FASTQ
2. run **Cutadapt**
3. run **FastQC** again on trimmed FASTQ
4. run **STAR 2-pass**
5. generate **samtools** summary statistics
6. generate **gene-level integer counts** with **featureCounts**
7. generate **gene-level expected counts** and **TPM** with **RSEM**
8. run **MultiQC once at the very end**

MultiQC is executed only after all sample-level steps complete successfully.

---

## 9. Main outputs

After successful completion, the main outputs are expected to be:

* `counts/gene_integer_counts.tsv`
* `abundance/gene_expected_counts.tsv`
* `abundance/gene_tpm.tsv`
* `multiqc/multiqc_report.html`

Additional outputs will include:

* trimmed FASTQ files
* BAM and BAI files
* FastQC reports
* STAR logs
* samtools stats
* per-sample logs
* metadata and status files

---

## 10. How to tell whether the run finished successfully

WulfRNA uses explicit status files so that completion is easy to judge.

### During execution

```text
status/RUNNING
```

### On success

```text
status/SUCCESS
status/finished_at.txt
status/summary.txt
```

### On failure

```text
status/FAILED
status/failed_step.txt
status/failed_sample.txt
status/summary.txt
```

### Quick check

To inspect pipeline state:

```bash
ls status/
cat status/summary.txt
```

If `SUCCESS` exists, the run completed successfully.
If `FAILED` exists, the run did not complete.

---

## 11. Recommended pre-run checklist

Before every run, check the following.

### A. Input files exist

```bash
ls fastq/*_R1.fastq.gz
ls fastq/*_R2.fastq.gz
```

### B. Required commands are on `PATH`

```bash
which fastqc
which cutadapt
which STAR
which samtools
which featureCounts
which rsem-calculate-expression
which multiqc
```

### C. Reference directory exists and is populated

```bash
ls /path/to/reference_dir
```

### D. Enough disk space is available

Large RNA-seq runs can consume substantial disk space for:

* trimmed FASTQ
* STAR temporary files
* BAM files
* RSEM intermediate files
* MultiQC artifacts

### E. Strandedness is known

Use the correct value for:

* `none`
* `forward`
* `reverse`

Do **not** guess silently. If strandedness is uncertain, verify it before running.

---

## 12. Example full session on bitwulf

This is only an example. Adjust module names and paths to your actual environment.

```bash
cd /data/fuki/projects/WulfRNA_test

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

fastqc --version
cutadapt --version
STAR --version
samtools --version
featureCounts -v
rsem-calculate-expression --version
multiqc --version

python3 -m venv .venv
source .venv/bin/activate
pip install -e /path/to/WulfRNA

wulfrna run /data/fuki/projects/WulfRNA_test \
  --reference /data/fuki/references/mm10_M21 \
  --stranded reverse \
  --threads 16
```

After completion:

```bash
cat status/summary.txt
ls multiqc/
ls counts/
ls abundance/
```

---

## 13. Common failure modes

### Missing R2 mate

Example:

```text
sampleA_R1.fastq.gz
```

but no matching `sampleA_R2.fastq.gz`

Result: the pipeline should fail before processing starts.

### Tool missing from `PATH`

Example:

```bash
which STAR
```

returns nothing.

Result: the pipeline is not runnable until the tool is installed or the correct module is loaded.

### Broken or incomplete reference directory

Result: the pipeline should fail before any sample-level work starts.

### Wrong strandedness

Result: counts and abundance may be biologically misleading even if the run technically succeeds.

---

## 14. Design philosophy

WulfRNA is intentionally narrow:

* paired-end bulk RNA-seq only
* single-server execution
* explicit status markers
* counts for SARTools / DESeq2
* TPM and expected counts for convenience
* final MultiQC report only after complete success

This is a lab pipeline, not a general workflow platform.

---

## 15. Development note

When implementing or extending WulfRNA, keep the following principles:

* do **not** hardcode server-specific executable paths
* rely on required tools being available in `PATH`
* fail early if dependencies are missing
* make completion obvious through `status/`
* keep the run model simple and auditable

---
