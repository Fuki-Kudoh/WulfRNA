import csv
from pathlib import Path

import pytest

from wulfrna.pipeline import GeneAnnotation, PipelineError, Sample, aggregate_star_counts, aggregate_transcript_quant


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


def write_star_file(workdir, sample, rows):
    source = workdir / "align" / "star" / sample / "ReadsPerGene.out.tab"
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_text(
        "N_unmapped\t0\t0\t0\nN_multimapping\t0\t0\t0\nN_noFeature\t0\t0\t0\nN_ambiguous\t0\t0\t0\n"
        + "".join(rows),
        encoding="utf-8",
    )
    return source


def run_star_aggregation(tmp_path, sample_rows, annotations, stranded="none"):
    (tmp_path / "abundance").mkdir()
    for sample, rows in sample_rows.items():
        write_star_file(tmp_path, sample, rows)
    aggregate_star_counts(
        tmp_path,
        [Sample(sample, Path("unused")) for sample in sample_rows],
        stranded,
        annotations,
    )


@pytest.mark.parametrize(
    "sample_rows, annotations, message",
    [
        (
            {"one": ["gene1\t1\t1\t1\n"], "two": ["gene2\t1\t1\t1\n"]},
            {"gene1": GeneAnnotation("G1", 100), "gene2": GeneAnnotation("G2", 100)},
            "gene sets are inconsistent",
        ),
        (
            {"one": ["gene1\t1\t1\t1\n", "gene1\t2\t2\t2\n"]},
            {"gene1": GeneAnnotation("G1", 100)},
            "Duplicate STAR gene_id",
        ),
        (
            {"one": ["gene1\tbad\t1\t1\n"]},
            {"gene1": GeneAnnotation("G1", 100)},
            "Non-integer STAR count",
        ),
        (
            {"one": ["gene1\t1\t1\n"]},
            {"gene1": GeneAnnotation("G1", 100)},
            "four columns",
        ),
        (
            {"one": ["gene1\t1\t1\t1\textra\n"]},
            {"gene1": GeneAnnotation("G1", 100)},
            "four columns",
        ),
        (
            {"one": ["unknown\t1\t1\t1\n"]},
            {"gene1": GeneAnnotation("G1", 100)},
            "absent from combined_gene_annotation",
        ),
        (
            {"one": ["gene1\t0\t0\t0\n"]},
            {"gene1": GeneAnnotation("G1", 100)},
            "Total STAR RPK is zero",
        ),
    ],
)
def test_star_input_validation(tmp_path, sample_rows, annotations, message):
    with pytest.raises(PipelineError, match=message) as error:
        run_star_aggregation(tmp_path, sample_rows, annotations)
    assert error.value.step == "aggregate"


def test_star_tpm_for_multiple_genes_and_samples(tmp_path):
    annotations = {
        "gene1": GeneAnnotation("G1", 1000),
        "gene2": GeneAnnotation("G2", 2000),
        "gene3": GeneAnnotation("G3", 500),
    }
    run_star_aggregation(
        tmp_path,
        {
            "sample_b": ["gene1\t0\t0\t0\n", "gene2\t4\t4\t4\n", "gene3\t1\t1\t1\n"],
            "sample_a": ["gene1\t2\t2\t2\n", "gene2\t2\t2\t2\n", "gene3\t0\t0\t0\n"],
        },
        annotations,
    )
    rows = list(csv.DictReader((tmp_path / "abundance/star_gene_tpm.tsv").open(), delimiter="\t"))
    assert [row["gene_id"] for row in rows] == ["gene1", "gene2", "gene3"]
    assert [row["sample_a"] for row in rows] == ["666666.666667", "333333.333333", "0.000000"]
    assert [row["sample_b"] for row in rows] == ["0.000000", "500000.000000", "500000.000000"]
    for sample in ["sample_a", "sample_b"]:
        assert abs(sum(float(row[sample]) for row in rows) - 1_000_000) <= 0.01


def test_star_tpm_validates_cumulative_serialized_rounding(tmp_path):
    gene_count = 1001
    annotations = {
        f"gene{i}": GeneAnnotation(f"G{i}", 101 + i % 17)
        for i in range(gene_count)
    }
    rows = [f"gene{i}\t{i % 11 + 1}\t0\t0\n" for i in range(gene_count)]
    run_star_aggregation(tmp_path, {"sample": rows}, annotations)

    with (tmp_path / "abundance/star_gene_tpm.tsv").open() as output:
        serialized = list(csv.DictReader(output, delimiter="\t"))
    serialized_sum = sum(float(row["sample"]) for row in serialized)
    assert abs(serialized_sum - 1_000_000) <= 0.01
