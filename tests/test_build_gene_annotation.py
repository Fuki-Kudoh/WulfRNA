import importlib.util
import sys
from pathlib import Path

import pytest


SCRIPT = Path(__file__).parents[1] / "scripts" / "build_gene_annotation.py"
SPEC = importlib.util.spec_from_file_location("build_gene_annotation", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
builder = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = builder
SPEC.loader.exec_module(builder)


def test_builds_sorted_exon_union_lengths_and_missing_names(tmp_path):
    gtf = tmp_path / "combined.gtf"
    gtf.write_text(
        'chr1\tsource\texon\t100\t200\t.\t+\t.\tgene_id "gene.2";\n'
        'chr1\tsource\texon\t1\t10\t.\t+\t.\tgene_id "gene.1"; gene_name "One";\n'
        'chr1\tsource\texon\t8\t20\t.\t+\t.\tgene_id "gene.1"; gene_name "One";\n'
        'chr1\tsource\texon\t30\t30\t.\t+\t.\tgene_id "gene.1"; gene_name "One";\n'
        'chr1\tsource\ttranscript\t1\t30\t.\t+\t.\tgene_id "ignored";\n',
        encoding="utf-8",
    )
    output = tmp_path / "combined_gene_annotation.tsv"

    builder.build_gene_annotation(gtf, output)

    assert output.read_text(encoding="utf-8").splitlines() == [
        "gene_id\tGeneName\tgene_length_bp",
        "gene.1\tOne\t21",
        "gene.2\tNA\t101",
    ]


def test_supports_realistic_gencode_unquoted_numeric_attributes(tmp_path):
    gtf = tmp_path / "gencode.gtf"
    gtf.write_text(
        'chr1\tHAVANA\texon\t100\t200\t.\t+\t.\t'
        'gene_id "ENSMUSG00000000001.1"; transcript_id "ENSMUST00000000001.1"; '
        'gene_name "Gnai3"; exon_number 1; level 2; tag "basic";\n'
        'chr1\tHAVANA\texon\t180\t250\t.\t+\t.\t'
        'gene_id "ENSMUSG00000000001.1"; transcript_id "ENSMUST00000000002.1"; '
        'gene_name "Gnai3"; exon_number 2; level 2; tag "basic";\n',
        encoding="utf-8",
    )
    output = tmp_path / "combined_gene_annotation.tsv"

    builder.build_gene_annotation(gtf, output)

    assert output.read_text(encoding="utf-8").splitlines() == [
        "gene_id\tGeneName\tgene_length_bp",
        "ENSMUSG00000000001.1\tGnai3\t151",
    ]


def test_rejects_gene_on_multiple_chromosomes(tmp_path):
    gtf = tmp_path / "combined.gtf"
    gtf.write_text(
        'chr1\ts\texon\t1\t10\t.\t+\t.\tgene_id "gene1";\n'
        'chr2\ts\texon\t20\t30\t.\t+\t.\tgene_id "gene1";\n',
        encoding="utf-8",
    )

    with pytest.raises(builder.AnnotationBuildError, match="multiple chromosomes"):
        builder.build_gene_annotation(gtf, tmp_path / "out.tsv")


@pytest.mark.parametrize(
    "line, message",
    [
        ("chr1\ts\texon\t1\t10\t.\t+\t.\n", "nine"),
        ('chr1\ts\texon\tbad\t10\t.\t+\t.\tgene_id "g";\n', "Non-integer"),
        ('chr1\ts\texon\t10\t1\t.\t+\t.\tgene_id "g";\n', "Invalid exon interval"),
        ('chr1\ts\texon\t1\t10\t.\t+\t.\tgene_name "G";\n', "no gene_id"),
    ],
)
def test_rejects_malformed_exons(tmp_path, line, message):
    gtf = tmp_path / "combined.gtf"
    gtf.write_text(line, encoding="utf-8")
    with pytest.raises(builder.AnnotationBuildError, match=message):
        builder.build_gene_annotation(gtf, tmp_path / "out.tsv")
