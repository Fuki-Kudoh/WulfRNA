# WulfRNA MVP Spec

## Scope
WulfRNA is a lightweight single-server bulk RNA-seq pipeline exposed as a packaged CLI.
Paired-end FASTQ input is the default, and single-end FASTQ input is supported with `--single-end` / `--SE`.

## CLI
Primary command:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded {none|forward|reverse} --threads N [--quantifier {salmon|kallisto}] [--aligner {none|star}] [--single-end|--SE] [--fragment-length FLOAT] [--fragment-sd FLOAT] [--dry-run] [--genome NAME]
```

- `WORKDIR` must contain `fastq/` using the selected input layout.
- `--reference` points to a reference root.
- `--genome` (optional) resolves references under `<reference>/<genome>/`.
- `--quantifier` defaults to `salmon` and selects the transcript quantification backend.
- `--aligner` defaults to `none`; `--aligner star` runs STAR from trimmed FASTQs and writes isolated outputs under `align/star/<sample>/` while preserving Salmon/kallisto quantification outputs.
- `--single-end` / `--SE` enables single-end input; paired-end is the default.
- `--fragment-length` and `--fragment-sd` are required only for single-end kallisto runs and must be positive.
- `--no-resume` disables automatic phase-level resume.
- `--force-from` forces rerun from one phase onward (`fastqc_raw`, `cutadapt`, `fastqc_trimmed`, `align`, `quant`, `aggregate`, `multiqc`).

## Required input layout
`WORKDIR/fastq/` must contain one of these layouts:

Paired-end default:
- `sample_id_R1.fastq.gz`
- `sample_id_R2.fastq.gz`

Single-end (`--single-end` / `--SE`):
- `sample_id.fastq.gz` or `sample_id_R1.fastq.gz`
- any `*_R2.fastq.gz` file is rejected in single-end mode.

## Required external tools
All must be in `PATH`:
- fastqc
- cutadapt
- multiqc
- salmon (if `--quantifier salmon`)
- kallisto (if `--quantifier kallisto`)
- STAR and samtools (if `--aligner star`)

## Required reference files (resolved reference directory)
Shared:
- `combined_tx2gene.tsv`
- `combined_gene_annotation.tsv` with unique, nonblank `gene_id`, `GeneName` (blank is normalized to `NA`), and positive integer `gene_length_bp`. Length is the exon-union length, so overlapping bases within a gene are counted once.

Salmon backend:
- `salmon_index/`

kallisto backend:
- `kallisto_index/combined_transcripts.kidx`

STAR aligner (`--aligner star`):
- `star_index/Genome`
- `star_index/SA`
- `star_index/SAindex`
- `star_index/genomeParameters.txt`
- `star_index/chrName.txt`
- `star_index/chrLength.txt`
- `star_index/chrNameLength.txt`

All required STAR index files must be non-empty.

### Reference preparation

Build the gene annotation from the same combined GTF used for the STAR index
and transcript-to-gene mapping:

```bash
python scripts/build_gene_annotation.py --gtf combined.gtf --output combined_gene_annotation.tsv
```

The builder considers exon features only, preserves GTF gene IDs, merges
overlapping/adjacent one-based intervals per gene, and emits their union length.
Missing `gene_name` becomes `NA`. Malformed exon records, invalid coordinates,
inconsistent gene names, and a gene occurring on multiple chromosomes are fatal.

## Pipeline behavior
1. Validate tools, references, and FASTQ inputs for the selected layout.
2. Record metadata (`samples.tsv`, run parameters, versions).
3. If `--dry-run`, stop after validation and metadata/status writing.
4. Otherwise run:
   - FastQC on raw FASTQ
   - Cutadapt trimming
   - FastQC on trimmed FASTQ
   - STAR alignment per sample when `--aligner star` is selected
   - Transcript quantification per sample (`salmon quant` or `kallisto quant`)
   - Aggregate transcript-level estimates to gene-level using `combined_tx2gene.tsv`
   - MultiQC over the work directory
5. Verify required final outputs exist.

### Phase-level checkpointing and resume

Completed phases write:
- `status/steps/fastqc_raw.done`
- `status/steps/cutadapt.done`
- `status/steps/fastqc_trimmed.done`
- `status/steps/align.done` (only when `--aligner star`)
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
- Includes `workdir`, `reference_dir`, `quantifier`, `aligner`, `stranded`, `layout`, `sample_ids`, reference files used, tx2gene and gene-annotation fingerprints, STAR index fingerprint when applicable, `created_at`, `updated_at`.

Manifest compatibility rules:
- sample set change => hard error (no silent resume),
- missing `layout` in old manifests is treated as `paired_end`,
- layout mismatch between existing and current manifest => hard error (prevents silent paired-end/single-end mixing),
- quantifier/reference_dir/stranded change => rerun quant and downstream,
- aligner or STAR index path/content change => rerun align and downstream when STAR alignment is enabled,
- `combined_tx2gene.tsv` fingerprint change => reuse quant, rerun aggregate + multiqc,
- `combined_gene_annotation.tsv` fingerprint change => reuse quant, rerun aggregate + multiqc,
- threads-only changes are resumable.

Failure handling:
- failed phase does not get a `.done` marker,
- `status/FAILED` and `status/failed_step.txt` are written,
- `status/failed_sample.txt` is written when sample context exists,
- previously completed phase markers remain for future resume.

## Output contracts
Primary outputs:
- `abundance/<quantifier>_gene_expected_counts.tsv` and `abundance/<quantifier>_gene_tpm.tsv`, whose leading columns are `gene_id` and `GeneName`
- legacy `abundance/gene_expected_counts.tsv` and `abundance/gene_tpm.tsv`, retaining the v0.2.x `gene_id` + sample-column format
- `multiqc/multiqc_report.html`
- `logs/tx2gene_mapping_stats.tsv`

STAR outputs when `--aligner star`:
- `align/star/<sample>/Aligned.sortedByCoord.out.bam`
- `align/star/<sample>/Aligned.sortedByCoord.out.bam.bai`
- `align/star/<sample>/SJ.out.tab`
- `align/star/<sample>/ReadsPerGene.out.tab`
- `align/star/<sample>/Log.final.out`
- copied native file `abundance/star_gene_counts/<sample>.star.ReadsPerGene.out.tab`
- integer matrix `abundance/star_gene_counts.tsv`
- `abundance/star_gene_tpm.tsv`, WulfRNA-calculated gene-level TPM from the strandedness-selected STAR count column and exon-union gene lengths

STAR summary rows are excluded from matrices but retained in copied native files. STAR TPM is `count / length_kb`, normalized per sample to one million; it is distinct from Salmon/kallisto transcript-quantifier TPM. Values are written with six decimal places, and the sum of those serialized values must be within an absolute tolerance of `0.01` of one million for every sample.

Status markers (`WORKDIR/status/`):
- `RUNNING` while active
- `SUCCESS` on complete run success
- `DRY_RUN_OK` on dry-run success
- `FAILED` on failure (+ `failed_step.txt`, optional `failed_sample.txt`)
- `summary.txt` and `finished_at.txt` on terminal states
- `steps/*.done` per completed phase
- `manifest.json` for resume compatibility
