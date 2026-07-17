import pytest

from wulfrna.pipeline import (
    PipelineError,
    STAR_INDEX_REQUIRED_FILES,
    compare_manifest,
    fingerprint_star_index,
    validate_reference,
    validate_star_index,
)


def make_salmon_reference(tmp_path):
    (tmp_path / "combined_tx2gene.tsv").write_text("tx1\tgene1\n", encoding="utf-8")
    (tmp_path / "salmon_index").mkdir()
    return tmp_path


def make_kallisto_reference(tmp_path):
    (tmp_path / "combined_tx2gene.tsv").write_text("tx1\tgene1\n", encoding="utf-8")
    kallisto_dir = tmp_path / "kallisto_index"
    kallisto_dir.mkdir()
    (kallisto_dir / "combined_transcripts.kidx").write_text("index\n", encoding="utf-8")
    return tmp_path


def make_complete_star_index(reference_dir):
    star_index = reference_dir / "star_index"
    star_index.mkdir()
    for name in STAR_INDEX_REQUIRED_FILES:
        (star_index / name).write_text(f"{name}\n", encoding="utf-8")
    return star_index


def test_validate_reference_default_aligner_does_not_require_star_index(tmp_path):
    ref = make_salmon_reference(tmp_path)

    refs = validate_reference(ref, "salmon")

    assert refs["tx2gene"] == ref / "combined_tx2gene.tsv"
    assert refs["index"] == ref / "salmon_index"
    assert "star_index" not in refs


def test_validate_star_index_missing_directory(tmp_path):
    star_index = tmp_path / "star_index"

    with pytest.raises(PipelineError) as exc:
        validate_star_index(star_index)

    assert exc.value.step == "reference_check"
    assert str(exc.value) == f"STAR index directory not found: {star_index}"


def test_validate_star_index_empty_directory(tmp_path):
    star_index = tmp_path / "star_index"
    star_index.mkdir()

    with pytest.raises(PipelineError) as exc:
        validate_star_index(star_index)

    assert exc.value.step == "reference_check"
    assert "STAR index is missing required non-empty files" in str(exc.value)
    for name in STAR_INDEX_REQUIRED_FILES:
        assert str(star_index / name) in str(exc.value)


def test_validate_star_index_one_missing_required_file(tmp_path):
    star_index = make_complete_star_index(tmp_path)
    missing_file = star_index / "SAindex"
    missing_file.unlink()

    with pytest.raises(PipelineError) as exc:
        validate_star_index(star_index)

    assert exc.value.step == "reference_check"
    assert str(missing_file) in str(exc.value)


def test_validate_star_index_one_zero_byte_required_file(tmp_path):
    star_index = make_complete_star_index(tmp_path)
    empty_file = star_index / "Genome"
    empty_file.write_bytes(b"")

    with pytest.raises(PipelineError) as exc:
        validate_star_index(star_index)

    assert exc.value.step == "reference_check"
    assert str(empty_file) in str(exc.value)


def test_validate_star_index_complete_index(tmp_path):
    star_index = make_complete_star_index(tmp_path)

    assert validate_star_index(star_index) == star_index


def test_validate_reference_salmon_star(tmp_path):
    ref = make_salmon_reference(tmp_path)
    star_index = make_complete_star_index(ref)

    refs = validate_reference(ref, "salmon", "star")

    assert refs["index"] == ref / "salmon_index"
    assert refs["star_index"] == star_index


def test_validate_reference_kallisto_star(tmp_path):
    ref = make_kallisto_reference(tmp_path)
    star_index = make_complete_star_index(ref)

    refs = validate_reference(ref, "kallisto", "star")

    assert refs["index"] == ref / "kallisto_index" / "combined_transcripts.kidx"
    assert refs["star_index"] == star_index


def test_validate_reference_rejects_unknown_aligner(tmp_path):
    ref = make_salmon_reference(tmp_path)

    with pytest.raises(PipelineError) as exc:
        validate_reference(ref, "salmon", "hisat2")

    assert exc.value.step == "reference_check"
    assert str(exc.value) == "Unsupported aligner: hisat2"


def test_compare_manifest_forces_align_when_star_index_content_changes(tmp_path):
    ref = make_salmon_reference(tmp_path)
    star_index = make_complete_star_index(ref)
    old_fingerprint = fingerprint_star_index(star_index)
    (star_index / "chrName.txt").write_text("changed chr names\n", encoding="utf-8")
    new_fingerprint = fingerprint_star_index(star_index)

    base_manifest = {
        "sample_ids": ["sample"],
        "layout": "paired_end",
        "quantifier": "salmon",
        "reference_dir": str(ref),
        "stranded": "reverse",
        "aligner": "star",
        "reference_files": {"star_index": str(star_index)},
        "tx2gene_fingerprint": {"sha256": "same"},
    }
    existing = {**base_manifest, "star_index_fingerprint": old_fingerprint}
    current = {**base_manifest, "star_index_fingerprint": new_fingerprint}

    forced_phase, reasons = compare_manifest(existing, current)

    assert forced_phase == "align"
    assert "star_index fingerprint changed" in reasons


def test_fingerprint_star_index_does_not_hash_large_index_file_contents(monkeypatch, tmp_path):
    ref = make_salmon_reference(tmp_path)
    star_index = make_complete_star_index(ref)
    original_open = type(star_index).open
    large_names = {"Genome", "SA", "SAindex"}

    def guarded_open(self, *args, **kwargs):
        if self.name in large_names:
            raise AssertionError(f"large STAR index file was opened for hashing: {self.name}")
        return original_open(self, *args, **kwargs)

    monkeypatch.setattr(type(star_index), "open", guarded_open)

    fingerprint = fingerprint_star_index(star_index)

    for name in large_names:
        assert "sha256" not in fingerprint["files"][name]
        assert set(fingerprint["files"][name]) == {"path", "size", "mtime_ns"}
    for name in {"genomeParameters.txt", "chrName.txt", "chrLength.txt", "chrNameLength.txt"}:
        assert "sha256" in fingerprint["files"][name]
