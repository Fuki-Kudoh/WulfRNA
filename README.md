# WulfRNA
<p align="center">
  <img src="assets/wulfrna_logo.png" alt="WulfRNA logo" width="500">
</p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21409117.svg)](https://doi.org/10.5281/zenodo.21409117)

WulfRNA is a packaged CLI for a focused bulk RNA-seq workflow on paired-end (default) or single-end FASTQ input.

Current outputs:
- source-qualified, gene-level expected counts and TPM from Salmon or kallisto
- optional STAR count and WulfRNA-calculated gene-level TPM matrices
- optional coordinate-sorted STAR BAM and BAI
- STAR splice-junction and gene-count outputs
- final MultiQC report

Pipeline steps: FastQC (raw) → Cutadapt → FastQC (trimmed) →
optional STAR alignment → transcript quantification (Salmon or kallisto) →
gene-level aggregation → MultiQC.

Public-readiness note: WulfRNA v0.2.2 is a lightweight single-server
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
pip install .
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
- shared: `combined_gene_annotation.tsv`, with columns `gene_id`, `GeneName`, and positive integer `gene_length_bp` (the union length of annotated exons, counting overlaps once)
- Salmon backend: `salmon_index/`
- kallisto backend: `kallisto_index/combined_transcripts.kidx`
- STAR aligner (when `--aligner star`): `star_index/` containing non-empty `Genome`, `SA`, `SAindex`, `genomeParameters.txt`, `chrName.txt`, `chrLength.txt`, and `chrNameLength.txt`

### Build the shared gene annotation

Generate `combined_gene_annotation.tsv` from the same GTF used to build the
STAR index and `combined_tx2gene.tsv`:

```bash
python scripts/build_gene_annotation.py \
  --gtf /path/to/combined.gtf \
  --output /path/to/reference/combined_gene_annotation.tsv
```

The builder reads exon records, preserves each GTF `gene_id` exactly, merges
overlapping or adjacent exon intervals within a gene, and writes their union
length as `gene_length_bp`. It uses `gene_name` for `GeneName` and writes `NA`
when no symbol exists. A gene assigned to multiple chromosomes, inconsistent
gene names, malformed records, or invalid coordinates cause a clear failure.
Use the identical combined GTF for all reference artifacts so gene IDs match
between this table, `combined_tx2gene.tsv`, and STAR.

## 4) Run command

```bash
wulfrna run WORKDIR --reference /path/to/reference --stranded reverse --threads 16 \
  --adapter-r1 R1_SEQUENCE --adapter-r2 R2_SEQUENCE
```

`R1_SEQUENCE` and `R2_SEQUENCE` are placeholders. Replace them with the regular
3′ adapter sequences specified by the library preparation protocol.

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
- `--adapter-r1 SEQUENCE` (repeatable): regular 3′ adapter candidate for read 1
- `--adapter-r2 SEQUENCE` (repeatable, paired-end only): regular 3′ adapter candidate for read 2
- `--no-adapter-trimming`: explicitly preserve the Cutadapt phase without removing adapters
- `--fragment-length FLOAT` and `--fragment-sd FLOAT` (single-end kallisto only; both required and must be `> 0`)
- `--min-mapping-rate FLOAT` (optional, default `0.90`): minimum acceptable tx2gene transcript mapping rate per sample (`0.0-1.0`)
- `--dry-run` (optional): validate tools/reference/inputs, write separate dry-run metadata/status, then exit
- `--genome NAME` (optional): resolve references as `<reference>/<NAME>/...`
- `--no-resume` (optional): disable phase-level resume and rerun all phases
- `--force-from {fastqc_raw,cutadapt,fastqc_trimmed,align,quant,aggregate,multiqc}` (optional): force rerun from the selected phase onward

Backward-compatibility note: legacy invocation without explicit `run` is still accepted, but `wulfrna run ...` is the intended interface.

### Adapter configuration (required)

Every run must choose exactly one explicit trimming strategy:

- Paired-end: provide at least one `--adapter-r1` and at least one `--adapter-r2`.
- Single-end: provide at least one `--adapter-r1`; `--adapter-r2` is rejected.
- Intentional adapter-free processing: provide `--no-adapter-trimming` and no adapter sequences.

Repeat an adapter option to test multiple protocol-appropriate candidates:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 8 \
  --adapter-r1 R1_SEQUENCE_1 --adapter-r1 R1_SEQUENCE_2 \
  --adapter-r2 R2_SEQUENCE_1 --adapter-r2 R2_SEQUENCE_2
```

Single-end example:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 8 \
  --single-end --adapter-r1 R1_SEQUENCE
```

Sequences are case-normalized and must contain only supported IUPAC DNA codes.
Advanced Cutadapt expressions, linked/anchored adapters, and presets are not
supported. Choose sequences from the library preparation protocol: a FastQC
adapter label does not establish the correct complete sequence, and WulfRNA does
not silently substitute a generic Illumina or Nextera adapter.

**Migration requirement:** runs that previously omitted adapter options must now
provide explicit sequences or `--no-adapter-trimming`. Existing work directories
whose manifests do not contain verified trimming metadata are invalidated from
Cutadapt onward; raw FastQC output remains reusable.

## 5) Minimal smoke test plan

```bash
pip install -e .
wulfrna --help
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --no-adapter-trimming --dry-run
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --quantifier kallisto --no-adapter-trimming --dry-run
wulfrna run <workdir> --reference <reference_dir> --stranded reverse --threads 4 --aligner star --no-adapter-trimming --dry-run
```

Notes:
- These smoke commands use `--no-adapter-trimming` intentionally to validate setup only; production runs should use protocol-specific adapter sequences when removal is required.
- For dry-run to pass, required binaries must be in `PATH`, references must be complete for the selected quantifier and optional aligner, and `<workdir>/fastq` must contain valid FASTQ inputs for the selected layout (`paired_end` or `single_end`). A dry run reports resume invalidation without replacing the authoritative execution manifest, status, `logs/run_parameters.txt`, or `metadata/versions.txt`; its invocation and versions are written to `logs/dry_run_parameters.txt` and `metadata/dry_run_versions.txt`, and its status is written under `status/dry_run/`.
- Single-end examples:
  - Salmon: `wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 4 --single-end --adapter-r1 R1_SEQUENCE --dry-run`
  - kallisto: `wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 4 --single-end --quantifier kallisto --fragment-length 200 --fragment-sd 20 --adapter-r1 R1_SEQUENCE --dry-run`

## 6) Main outputs and status markers

Expected primary outputs on full success:
- `abundance/<quantifier>_gene_expected_counts.tsv`
- `abundance/<quantifier>_gene_tpm.tsv`
- `abundance/gene_expected_counts.tsv` (legacy v0.2.x compatibility output without `GeneName`)
- `abundance/gene_tpm.tsv` (legacy v0.2.x compatibility output without `GeneName`)
- `multiqc/multiqc_report.html`
- `logs/tx2gene_mapping_stats.tsv`

Additional outputs when `--aligner star` is selected (per sample):
- `align/star/<sample>/Aligned.sortedByCoord.out.bam`
- `align/star/<sample>/Aligned.sortedByCoord.out.bam.bai`
- `align/star/<sample>/SJ.out.tab`
- `align/star/<sample>/ReadsPerGene.out.tab`
- `align/star/<sample>/Log.final.out`
- `abundance/star_gene_counts/<sample>.star.ReadsPerGene.out.tab` (metadata-preserving copy)
- `abundance/star_gene_counts.tsv`
- `abundance/star_gene_tpm.tsv`

New matrices begin with `gene_id`, `GeneName`, then deterministically sorted sample columns. Salmon/kallisto TPM is summed from transcript-quantifier output. In contrast, `star_gene_tpm.tsv` is **WulfRNA-calculated gene-level TPM from STAR counts** using exon-union lengths from `combined_gene_annotation.tsv`; it is not a native STAR output. STAR's count column is selected from `--stranded` (`none`: column 2, `forward`: column 3, `reverse`: column 4). STAR TPM values are serialized with six decimal places, then each serialized sample column is required to sum to 1,000,000 within an absolute tolerance of `0.01`.

Matrix behavior note:
- Gene-level matrices include only genes observed in the transcript quantification input (unobserved zero-only genes are not emitted).

Status files in `WORKDIR/status/`:
- during run: `RUNNING`
- full success: `SUCCESS`, `finished_at.txt`, `summary.txt`
- dry-run state: `dry_run/RUNNING`, then `dry_run/DRY_RUN_OK` or `dry_run/FAILED`, with dry-run-specific details in the same subdirectory
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
- If the adapter mode or either ordered adapter list changed, if a legacy manifest lacks trimming metadata, or if completed Cutadapt/downstream checkpoints have no manifest: `cutadapt` and every downstream phase rerun; `fastqc_raw` remains reusable.
- If `quantifier`, `reference_dir`, or `stranded` changed: `quant` and downstream phases rerun.
- If `aligner` or STAR index content changed: `align` and downstream phases rerun when STAR alignment is enabled.
- If `combined_tx2gene.tsv` fingerprint changed: `aggregate` and `multiqc` rerun.
- If `combined_gene_annotation.tsv` fingerprint changed: `aggregate` and `multiqc` rerun.
- If only thread count changed: resume is allowed.

To fully disable resume:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 16 \
  --adapter-r1 R1_SEQUENCE --adapter-r2 R2_SEQUENCE --no-resume
```

To explicitly rebuild trimmed reads and all dependent outputs with the same
configuration, use `--force-from cutadapt`. WulfRNA clears the Cutadapt and
downstream phase markers before recording new execution provenance, so an
interrupted rerun resumes at the first incomplete phase without accepting stale
outputs.

To rerun only aggregation/reporting after fixing `combined_tx2gene.tsv`:
- preferred: `--force-from aggregate`, or
- remove aggregate + multiqc outputs and rerun.

Example:

```bash
wulfrna run WORKDIR --reference REFDIR --stranded reverse --threads 16 \
  --adapter-r1 R1_SEQUENCE --adapter-r2 R2_SEQUENCE --force-from aggregate
```

## 8) Citation

If you use WulfRNA in your research, please cite:

> Kudoh, F. (2026). WulfRNA (Version 0.2.2) [Computer software].
> Zenodo. https://doi.org/10.5281/zenodo.21409117

Citation metadata is also available in `CITATION.cff`.
