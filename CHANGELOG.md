# Changelog

## [0.2.4] - 2026-10-21

### Fixed

- Cutadapt now removes reads shorter than 1 nt after adapter processing, preventing zero-length reads from reaching STAR.
- Paired-end filtering explicitly uses `--pair-filter=any`, so both mates are discarded when either mate is empty and retained pairs remain synchronized.

### Changed

- Run manifests now record the effective read-filtering policy. Manifests from v0.2.3 and earlier that lack this metadata automatically invalidate Cutadapt and all downstream phases while preserving raw FastQC.
- `--no-adapter-trimming` continues to disable adapter removal only; nonempty-read safety filtering remains active.

### Validation

- Added real-Cutadapt paired-end and single-end tests covering pre-existing empty reads and reads made empty by adapter removal.
- Added resume and dry-run regression coverage for legacy manifests without read-filtering metadata.

## [0.2.2] - 2026-07-21

### Added

- Source-qualified Salmon and kallisto matrices with gene symbols.
- Validated and fingerprinted `combined_gene_annotation.tsv` reference metadata.
- Collected STAR gene-count files, an integer gene-count matrix, and WulfRNA-calculated gene-level TPM based on exon-union lengths.
- Added a reproducible GTF exon-union builder for `combined_gene_annotation.tsv`.

### Changed

- Aggregate resume and final validation now require all applicable canonical outputs.
- Legacy unqualified matrices remain available for v0.2.x compatibility.
- STAR TPM validation now checks the serialized six-decimal values, and resume/error-path coverage includes v0.2.1 migration.

## [0.2.1] - 2026-07-17

### Added

- Added `CITATION.cff` with software citation metadata for GitHub and
  Zenodo integration.

### Changed

- Updated package version to 0.2.1 for the citation-enabled archived
  release.

## [0.2.0] - 2026-07-17

### Added

- Optional STAR alignment with `--aligner star`.
- Paired-end and single-end STAR support.
- Two-pass STAR alignment.
- Coordinate-sorted BAM and BAM index generation.
- Per-sample `SJ.out.tab`, `ReadsPerGene.out.tab`, and `Log.final.out`.
- STAR index validation and lightweight fingerprinting.
- Alignment-aware resume and `--force-from align`.
- STAR and samtools version capture.
- Automated pytest and fake-tool integration tests.

### Changed

- STAR outputs are isolated under `align/star/<sample>/`.
- Salmon and kallisto abundance outputs remain unchanged.

### Validation

- Validated using a 500,000-pair subset of GSE157878 BMDM RNA-seq
  against mm10.
- BAM passed `samtools quickcheck`.
- Unique mapping rate was 68.51%.
- A repeated invocation resumed successfully and skipped all phases.
