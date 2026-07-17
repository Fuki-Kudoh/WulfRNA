import pytest

from wulfrna.pipeline import PipelineError, validate_reference


def make_salmon_reference(tmp_path):
    (tmp_path / "combined_tx2gene.tsv").write_text("tx1\tgene1\n", encoding="utf-8")
    (tmp_path / "salmon_index").mkdir()
    return tmp_path


def test_validate_reference_default_aligner_does_not_require_star_index(tmp_path):
    ref = make_salmon_reference(tmp_path)

    refs = validate_reference(ref, "salmon")

    assert refs["tx2gene"] == ref / "combined_tx2gene.tsv"
    assert refs["index"] == ref / "salmon_index"
    assert "star_index" not in refs


def test_validate_reference_star_aligner_requires_star_index(tmp_path):
    ref = make_salmon_reference(tmp_path)

    with pytest.raises(PipelineError) as exc:
        validate_reference(ref, "salmon", "star")

    assert exc.value.step == "reference_check"
    assert "Reference directory is missing required files" in str(exc.value)
    assert "star_index" in str(exc.value)


def test_validate_reference_star_aligner_returns_star_index(tmp_path):
    ref = make_salmon_reference(tmp_path)
    (ref / "star_index").mkdir()

    refs = validate_reference(ref, "salmon", "star")

    assert refs["star_index"] == ref / "star_index"


def test_validate_reference_rejects_unknown_aligner(tmp_path):
    ref = make_salmon_reference(tmp_path)

    with pytest.raises(PipelineError) as exc:
        validate_reference(ref, "salmon", "hisat2")

    assert exc.value.step == "reference_check"
    assert str(exc.value) == "Unsupported aligner: hisat2"
