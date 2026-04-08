# WulfRNA MVP Spec

## Scope
WulfRNA is a single-server paired-end bulk RNA-seq pipeline exposed as a packaged CLI.

## CLI
Primary command:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded {none|forward|reverse} --threads N [--dry-run] [--genome NAME]
```

- `WORKDIR` must contain `fastq/` with `*_R1.fastq.gz` + matching `*_R2.fastq.gz`.
- `--reference` points to a reference root.
- `--genome` (optional) resolves references under `<reference>/<genome>/`.

## Required external tools
All must be in `PATH`:
- fastqc
- cutadapt
- STAR
- samtools
- featureCounts
- rsem-calculate-expression
- multiqc

## Required reference files (resolved reference directory)
- `star_index/`
- `annotation.gtf`
- `featurecounts.gtf`
- `rsem_reference.grp`
- `rsem_reference.ti`

## Pipeline behavior
1. Validate tools, references, and sample pairing.
2. Record metadata (`samples.tsv`, run parameters, versions).
3. If `--dry-run`, stop after validation and metadata/status writing.
4. Otherwise run:
   - FastQC on raw FASTQ
   - Cutadapt trimming
   - FastQC on trimmed FASTQ
   - STAR first pass per sample
   - STAR second pass per sample with combined splice junctions
   - samtools index/flagstat/stats
   - featureCounts per sample, then aggregate matrix
   - RSEM per sample, then aggregate expected-count and TPM matrices
   - MultiQC over the work directory
5. Verify required final outputs exist.

## Output contracts
Primary outputs:
- `counts/gene_integer_counts.tsv`
- `abundance/gene_expected_counts.tsv`
- `abundance/gene_tpm.tsv`
- `multiqc/multiqc_report.html`

Status markers (`WORKDIR/status/`):
- `RUNNING` while active
- `SUCCESS` on complete run success
- `DRY_RUN_OK` on dry-run success
- `FAILED` on failure (+ `failed_step.txt`, optional `failed_sample.txt`)
- `summary.txt` and `finished_at.txt` on terminal states
