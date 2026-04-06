# WulfRNA

Bulk RNA-seq pipeline MVP for bitwulf single-node execution.

## CLI

```bash
./rnaseq-run /path/to/workdir --reference /path/to/reference --genome hg38 --stranded reverse --threads 16
```

Required arguments:

- `workdir`
- `--reference`
- `--stranded` (`none|forward|reverse`)
- `--threads`

Optional arguments:

- `--genome`（`--reference/<genome>/` を参照する）

## Input layout

```text
workdir/
└── fastq/
    ├── <sample_id>_R1.fastq.gz
    └── <sample_id>_R2.fastq.gz
```

## Reference layout

`--genome` を指定しない場合、`--reference` 直下に次を配置します。指定した場合は `--reference/<genome>/` 配下に同じ構造を配置します。

- `star_index/`
- `annotation.gtf`
- `featurecounts.gtf`
- `rsem_reference` prefix with at least `rsem_reference.grp` and `rsem_reference.ti`

## Main outputs

- `counts/gene_integer_counts.tsv`
- `abundance/gene_expected_counts.tsv`
- `abundance/gene_tpm.tsv`
- `multiqc/multiqc_report.html`

The pipeline also writes `samples.tsv`, per-step logs under `logs/`, tool metadata under `metadata/versions.txt`, and sentinel files under `status/` (`RUNNING`, `SUCCESS`, `FAILED`).
