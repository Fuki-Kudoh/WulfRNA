# WulfRNA MVP Spec

## Scope
WulfRNA is a single-server paired-end bulk RNA-seq pipeline exposed as a packaged CLI.

## CLI
Primary command:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded {none|forward|reverse} --threads N [--quantifier {salmon|kallisto}] [--dry-run] [--genome NAME]
```

- `WORKDIR` must contain `fastq/` with `*_R1.fastq.gz` + matching `*_R2.fastq.gz`.
- `--reference` points to a reference root.
- `--genome` (optional) resolves references under `<reference>/<genome>/`.
- `--quantifier` defaults to `salmon` and selects transcript quant backend.
- `--no-resume` disables automatic phase-level resume.
- `--force-from` forces rerun from one phase onward (`fastqc_raw`, `cutadapt`, `fastqc_trimmed`, `quant`, `aggregate`, `multiqc`).

## Required external tools
All must be in `PATH`:
- fastqc
- cutadapt
- multiqc
- salmon (if `--quantifier salmon`)
- kallisto (if `--quantifier kallisto`)

## Required reference files (resolved reference directory)
Shared:
- `combined_tx2gene.tsv`

Salmon backend:
- `salmon_index/`

kallisto backend:
- `kallisto_index/combined_transcripts.kidx`

## Pipeline behavior
1. Validate tools, references, and sample pairing.
2. Record metadata (`samples.tsv`, run parameters, versions).
3. If `--dry-run`, stop after validation and metadata/status writing.
4. Otherwise run:
   - FastQC on raw FASTQ
   - Cutadapt trimming
   - FastQC on trimmed FASTQ
   - Transcript quantification per sample (`salmon quant` or `kallisto quant`)
   - Aggregate transcript-level estimates to gene-level using `combined_tx2gene.tsv`
   - MultiQC over the work directory
5. Verify required final outputs exist.

### Phase-level checkpointing and resume

Completed phases write:
- `status/steps/fastqc_raw.done`
- `status/steps/cutadapt.done`
- `status/steps/fastqc_trimmed.done`
- `status/steps/quant.done`
- `status/steps/aggregate.done`
- `status/steps/multiqc.done`

Each `.done` file contains an ISO timestamp.

Resume is enabled by default. A phase is skipped only if:
- done marker exists,
- expected outputs for that phase exist and are non-empty,
- run manifest compatibility allows reuse.

Machine-readable manifest:
- `status/manifest.json`
- Includes `workdir`, `reference_dir`, `quantifier`, `stranded`, `sample_ids`, reference files used, tx2gene fingerprint, `created_at`, `updated_at`.

Manifest compatibility rules:
- sample set change => hard error (no silent resume),
- quantifier/reference_dir/stranded change => rerun quant and downstream,
- `combined_tx2gene.tsv` fingerprint change => reuse quant, rerun aggregate + multiqc,
- threads-only changes are resumable.

Failure handling:
- failed phase does not get a `.done` marker,
- `status/FAILED` and `status/failed_step.txt` are written,
- `status/failed_sample.txt` is written when sample context exists,
- previously completed phase markers remain for future resume.

## Output contracts
Primary outputs:
- `abundance/gene_expected_counts.tsv`
- `abundance/gene_tpm.tsv`
- `multiqc/multiqc_report.html`

Status markers (`WORKDIR/status/`):
- `RUNNING` while active
- `SUCCESS` on complete run success
- `DRY_RUN_OK` on dry-run success
- `FAILED` on failure (+ `failed_step.txt`, optional `failed_sample.txt`)
- `summary.txt` and `finished_at.txt` on terminal states
- `steps/*.done` per completed phase
- `manifest.json` for resume compatibility
