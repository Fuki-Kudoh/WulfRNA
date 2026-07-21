import csv
from pathlib import Path

import pytest

from wulfrna.pipeline import GeneAnnotation, Sample, aggregate_star_counts, aggregate_transcript_quant


@pytest.mark.parametrize(
    "quantifier, header, row",
    [
        ("salmon", "Name\tNumReads\tTPM\n", "tx1\t2.5\t8.0\n"),
        ("kallisto", "target_id\test_counts\ttpm\n", "tx1\t2.5\t8.0\n"),
    ],
)
def test_source_qualified_and_legacy_transcript_matrices(tmp_path, quantifier, header, row):
    quant = tmp_path / "quant.tsv"
    quant.write_text(header + row, encoding="utf-8")
    canonical_counts = tmp_path / f"{quantifier}_gene_expected_counts.tsv"
    canonical_tpm = tmp_path / f"{quantifier}_gene_tpm.tsv"
    legacy_counts = tmp_path / "gene_expected_counts.tsv"
    legacy_tpm = tmp_path / "gene_tpm.tsv"
    aggregate_transcript_quant(
        {"sample": quant}, {"tx1": "gene1"}, quantifier,
        canonical_counts, canonical_tpm, tmp_path / "mapping.tsv", 0.9,
        {"gene1": GeneAnnotation("G1", 1000)}, legacy_counts, legacy_tpm,
    )
    assert canonical_counts.read_text().splitlines() == ["gene_id\tGeneName\tsample", "gene1\tG1\t2.500000"]
    assert canonical_tpm.read_text().splitlines() == ["gene_id\tGeneName\tsample", "gene1\tG1\t8.000000"]
    assert legacy_counts.read_text().splitlines() == ["gene_id\tsample", "gene1\t2.500000"]
    assert legacy_tpm.read_text().splitlines() == ["gene_id\tsample", "gene1\t8.000000"]


@pytest.mark.parametrize("stranded, expected", [("none", 2), ("forward", 3), ("reverse", 4)])
def test_star_strandedness_selects_expected_column(tmp_path, stranded, expected):
    source = tmp_path / "align" / "star" / "sample" / "ReadsPerGene.out.tab"
    source.parent.mkdir(parents=True)
    source.write_text(
        "N_unmapped\t0\t0\t0\nN_multimapping\t0\t0\t0\nN_noFeature\t0\t0\t0\nN_ambiguous\t0\t0\t0\n"
        "gene1\t2\t3\t4\n",
        encoding="utf-8",
    )
    (tmp_path / "abundance").mkdir()
    aggregate_star_counts(
        tmp_path, [Sample("sample", Path("unused"))], stranded,
        {"gene1": GeneAnnotation("G1", 1000)},
    )
    with (tmp_path / "abundance" / "star_gene_counts.tsv").open() as f:
        rows = list(csv.reader(f, delimiter="\t"))
    assert rows[1] == ["gene1", "G1", str(expected)]
    assert (tmp_path / "abundance" / "star_gene_tpm.tsv").read_text().splitlines()[1].endswith("1000000.000000")
