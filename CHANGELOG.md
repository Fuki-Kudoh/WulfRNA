# Changelog

## [0.2.2] - 2026-07-21

### Added

- Source-qualified Salmon and kallisto matrices with gene symbols.
- Validated and fingerprinted `combined_gene_annotation.tsv` reference metadata.
- Collected STAR gene-count files, an integer gene-count matrix, and WulfRNA-calculated gene-level TPM based on exon-union lengths.

### Changed

- Aggregate resume and final validation now require all applicable canonical outputs.
- Legacy unqualified matrices remain available for v0.2.x compatibility.

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
