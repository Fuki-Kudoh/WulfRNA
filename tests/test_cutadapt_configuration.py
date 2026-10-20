import gzip
import json
import shutil
import stat
import subprocess
from pathlib import Path

import pytest

from wulfrna.cli import main, parse_args
from wulfrna.pipeline import (
    LAYOUT_PAIRED,
    LAYOUT_SINGLE,
    PHASES,
    TRIMMING_ADAPTERS,
    TRIMMING_DISABLED,
    PipelineError,
    Sample,
    TrimmingConfig,
    build_cutadapt_command,
    compare_manifest,
    fingerprint_file,
    invalidate_phase_markers,
    load_manifest,
    mark_phase_done,
    should_run_phase,
    step_done_file,
    validate_cutadapt_report,
    write_manifest,
)


def base_args():
    return [
        "run",
        "work",
        "--reference",
        "ref",
        "--stranded",
        "reverse",
        "--threads",
        "4",
    ]


def test_explicit_paired_end_adapter_configuration_is_accepted_and_normalized():
    args = parse_args([*base_args(), "--adapter-r1", "acgtn", "--adapter-r2", "TTRYS"])

    assert args.adapter_r1 == ["ACGTN"]
    assert args.adapter_r2 == ["TTRYS"]
    assert args.no_adapter_trimming is False


def test_explicit_single_end_adapter_configuration_is_accepted():
    args = parse_args([*base_args(), "--single-end", "--adapter-r1", "ACGT"])

    assert args.adapter_r1 == ["ACGT"]
    assert args.adapter_r2 is None


def test_multiple_adapter_sequences_are_accepted():
    args = parse_args(
        [
            *base_args(),
            "--adapter-r1",
            "AAAA",
            "--adapter-r1",
            "CCCC",
            "--adapter-r2",
            "GGGG",
            "--adapter-r2",
            "TTTT",
        ]
    )

    assert args.adapter_r1 == ["AAAA", "CCCC"]
    assert args.adapter_r2 == ["GGGG", "TTTT"]


@pytest.mark.parametrize(
    "extra_args, message",
    [
        ([], "Adapter configuration is required"),
        (["--adapter-r1", "AAAA"], "--adapter-r1 and --adapter-r2"),
        (["--single-end", "--adapter-r1", "AAAA", "--adapter-r2", "CCCC"], "--adapter-r2 is not allowed"),
        (["--adapter-r1", "AAAA", "--adapter-r2", "CCCC", "--no-adapter-trimming"], "cannot be combined"),
        (["--adapter-r1", "ACGT;rm", "--adapter-r2", "CCCC"], "supported IUPAC DNA codes"),
        (["--adapter-r1", "", "--adapter-r2", "CCCC"], "supported IUPAC DNA codes"),
    ],
)
def test_invalid_adapter_configurations_are_rejected(extra_args, message, capsys):
    with pytest.raises(SystemExit) as exc:
        parse_args([*base_args(), *extra_args])

    assert exc.value.code == 2
    assert message in capsys.readouterr().err


def test_explicit_no_trimming_mode_is_accepted():
    args = parse_args([*base_args(), "--no-adapter-trimming"])

    assert args.no_adapter_trimming is True
    assert args.adapter_r1 is None
    assert args.adapter_r2 is None


def test_missing_adapter_configuration_is_rejected_before_workdir_state_is_created(tmp_path):
    workdir = tmp_path / "work"

    with pytest.raises(SystemExit):
        main(
            [
                "run",
                str(workdir),
                "--reference",
                str(tmp_path / "ref"),
                "--stranded",
                "reverse",
                "--threads",
                "2",
            ]
        )

    assert not workdir.exists()


def test_build_cutadapt_command_for_paired_end_preserves_paths_and_all_adapters(tmp_path):
    sample = Sample("sample", tmp_path / "fastq/sample_R1.fastq.gz", tmp_path / "fastq/sample_R2.fastq.gz")
    trimming = TrimmingConfig(TRIMMING_ADAPTERS, ("AAAA", "CCCC"), ("GGGG", "TTTT"))

    cmd = build_cutadapt_command(tmp_path, sample, 8, LAYOUT_PAIRED, trimming)

    assert [cmd[i + 1] for i, value in enumerate(cmd) if value == "-a"] == ["AAAA", "CCCC"]
    assert [cmd[i + 1] for i, value in enumerate(cmd) if value == "-A"] == ["GGGG", "TTTT"]
    assert cmd[cmd.index("-o") + 1] == str(tmp_path / "trimmed/sample_R1.trimmed.fastq.gz")
    assert cmd[cmd.index("-p") + 1] == str(tmp_path / "trimmed/sample_R2.trimmed.fastq.gz")
    assert cmd[cmd.index("--json") + 1] == str(tmp_path / "qc/cutadapt/sample.cutadapt.json")
    assert cmd[-2:] == [str(sample.r1), str(sample.r2)]


def test_build_cutadapt_command_for_single_end_uses_only_r1(tmp_path):
    sample = Sample("sample", tmp_path / "fastq/sample.fastq.gz")
    trimming = TrimmingConfig(TRIMMING_ADAPTERS, ("AAAA", "CCCC"), ())

    cmd = build_cutadapt_command(tmp_path, sample, 3, LAYOUT_SINGLE, trimming)

    assert [cmd[i + 1] for i, value in enumerate(cmd) if value == "-a"] == ["AAAA", "CCCC"]
    assert "-A" not in cmd
    assert "-p" not in cmd
    assert cmd[cmd.index("-o") + 1] == str(tmp_path / "trimmed/sample.trimmed.fastq.gz")
    assert cmd[-1] == str(sample.r1)


def test_build_cutadapt_command_for_disabled_trimming_omits_adapter_options(tmp_path):
    sample = Sample("sample", tmp_path / "fastq/sample_R1.fastq.gz", tmp_path / "fastq/sample_R2.fastq.gz")
    trimming = TrimmingConfig(TRIMMING_DISABLED, (), ())

    cmd = build_cutadapt_command(tmp_path, sample, 4, LAYOUT_PAIRED, trimming)

    assert "-a" not in cmd
    assert "-A" not in cmd


def write_report(path: Path, adapters_r1, adapters_r2):
    path.write_text(
        json.dumps({"adapters_read1": adapters_r1, "adapters_read2": adapters_r2}),
        encoding="utf-8",
    )


def adapter_definition(sequence: str, total_matches: int = 0):
    return {"name": "adapter", "sequence": sequence, "total_matches": total_matches}


def test_cutadapt_report_accepts_expected_adapters_with_zero_matches(tmp_path):
    report = tmp_path / "report.json"
    write_report(report, [adapter_definition("AAAA", 0)], [adapter_definition("CCCC", 0)])

    validate_cutadapt_report(
        report,
        TrimmingConfig(TRIMMING_ADAPTERS, ("AAAA",), ("CCCC",)),
        LAYOUT_PAIRED,
    )


def test_cutadapt_report_rejects_empty_definitions_when_trimming_requested(tmp_path):
    report = tmp_path / "report.json"
    write_report(report, [], [])

    with pytest.raises(PipelineError, match="R1 adapter definitions do not match"):
        validate_cutadapt_report(
            report,
            TrimmingConfig(TRIMMING_ADAPTERS, ("AAAA",), ("CCCC",)),
            LAYOUT_PAIRED,
        )


def test_cutadapt_report_accepts_empty_definitions_when_trimming_disabled(tmp_path):
    report = tmp_path / "report.json"
    write_report(report, [], [])

    validate_cutadapt_report(report, TrimmingConfig(TRIMMING_DISABLED, (), ()), LAYOUT_PAIRED)


@pytest.mark.parametrize("contents, message", [("not json", "is invalid"), ("{}", "is missing adapters_read1")])
def test_cutadapt_report_rejects_malformed_or_missing_data(tmp_path, contents, message):
    report = tmp_path / "report.json"
    report.write_text(contents, encoding="utf-8")

    with pytest.raises(PipelineError, match=message):
        validate_cutadapt_report(report, TrimmingConfig(TRIMMING_DISABLED, (), ()), LAYOUT_PAIRED)


def test_cutadapt_report_rejects_missing_file(tmp_path):
    with pytest.raises(PipelineError, match="missing or empty"):
        validate_cutadapt_report(
            tmp_path / "missing.json",
            TrimmingConfig(TRIMMING_DISABLED, (), ()),
            LAYOUT_PAIRED,
        )


def manifest(trimming):
    return {
        "sample_ids": ["sample"],
        "layout": LAYOUT_PAIRED,
        "quantifier": "salmon",
        "reference_dir": "/ref",
        "stranded": "reverse",
        "aligner": "none",
        "reference_files": {},
        "tx2gene_fingerprint": {"sha256": "tx"},
        "gene_annotation_fingerprint": {"sha256": "annotation"},
        "trimming": trimming,
    }


def trimming_manifest(r1=("AAAA",), r2=("CCCC",), mode=TRIMMING_ADAPTERS):
    return {"mode": mode, "adapters_r1": list(r1), "adapters_r2": list(r2)}


def test_unchanged_trimming_configuration_preserves_resume():
    existing = manifest(trimming_manifest())

    forced_phase, reasons = compare_manifest(existing, manifest(trimming_manifest()))

    assert forced_phase is None
    assert reasons == []


@pytest.mark.parametrize(
    "changed, reason",
    [
        (trimming_manifest(r1=("GGGG",)), "trimming configuration changed"),
        (trimming_manifest(r2=("TTTT",)), "trimming configuration changed"),
        (trimming_manifest(r1=(), r2=(), mode=TRIMMING_DISABLED), "trimming configuration changed"),
    ],
)
def test_changed_trimming_configuration_invalidates_from_cutadapt(changed, reason):
    forced_phase, reasons = compare_manifest(manifest(trimming_manifest()), manifest(changed))

    assert forced_phase == "cutadapt"
    assert reason in reasons


def test_legacy_manifest_without_trimming_metadata_invalidates_from_cutadapt():
    existing = manifest(trimming_manifest())
    existing.pop("trimming")

    forced_phase, reasons = compare_manifest(existing, manifest(trimming_manifest()))

    assert forced_phase == "cutadapt"
    assert "legacy manifest lacks verified trimming configuration" in reasons


def test_interrupted_configuration_change_cannot_reuse_downstream_markers(tmp_path):
    status_steps = tmp_path / "status/steps"
    status_steps.mkdir(parents=True)
    for phase in PHASES:
        mark_phase_done(tmp_path, phase)

    current = manifest(trimming_manifest(r1=("GGGG",)))
    invalidate_phase_markers(tmp_path, "cutadapt")
    write_manifest(tmp_path, current)

    assert step_done_file(tmp_path, "fastqc_raw").exists()
    assert all(not step_done_file(tmp_path, phase).exists() for phase in PHASES[1:])
    assert load_manifest(tmp_path)["trimming"] == current["trimming"]
    sample = Sample("sample", tmp_path / "fastq/sample_R1.fastq.gz", tmp_path / "fastq/sample_R2.fastq.gz")
    assert should_run_phase(tmp_path, "cutadapt", [sample], "salmon", LAYOUT_PAIRED, False, None)


def write_executable(path: Path):
    path.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def test_dry_run_with_changed_trimming_does_not_overwrite_manifest_or_markers(monkeypatch, tmp_path, capsys):
    from wulfrna.pipeline import execute

    workdir = tmp_path / "work"
    reference = tmp_path / "reference"
    bin_dir = tmp_path / "bin"
    (workdir / "fastq").mkdir(parents=True)
    (workdir / "fastq/sample_R1.fastq.gz").write_text("reads\n", encoding="utf-8")
    (workdir / "fastq/sample_R2.fastq.gz").write_text("reads\n", encoding="utf-8")
    reference.mkdir()
    (reference / "combined_tx2gene.tsv").write_text("tx1\tgene1\n", encoding="utf-8")
    (reference / "combined_gene_annotation.tsv").write_text(
        "gene_id\tGeneName\tgene_length_bp\ngene1\tGene1\t1000\n",
        encoding="utf-8",
    )
    (reference / "salmon_index").mkdir()
    bin_dir.mkdir()
    for tool in ["fastqc", "cutadapt", "multiqc", "salmon"]:
        write_executable(bin_dir / tool)
    monkeypatch.setenv("PATH", f"{bin_dir}:/usr/bin")

    existing = {
        **manifest(trimming_manifest()),
        "workdir": str(workdir.resolve()),
        "reference_dir": str(reference.resolve()),
        "reference_files": {
            "index": str(reference / "salmon_index"),
            "combined_tx2gene_tsv": str(reference / "combined_tx2gene.tsv"),
            "combined_gene_annotation_tsv": str(reference / "combined_gene_annotation.tsv"),
        },
        "tx2gene_fingerprint": fingerprint_file(reference / "combined_tx2gene.tsv"),
        "gene_annotation_fingerprint": fingerprint_file(reference / "combined_gene_annotation.tsv"),
    }
    (workdir / "status/steps").mkdir(parents=True)
    for phase in [phase for phase in PHASES if phase != "align"]:
        mark_phase_done(workdir, phase)
    write_manifest(workdir, existing)
    manifest_before = (workdir / "status/manifest.json").read_bytes()
    authoritative_params = b"completed run parameters\n"
    authoritative_versions = b"completed run versions\n"
    (workdir / "logs").mkdir(exist_ok=True)
    (workdir / "metadata").mkdir(exist_ok=True)
    (workdir / "logs/run_parameters.txt").write_bytes(authoritative_params)
    (workdir / "metadata/versions.txt").write_bytes(authoritative_versions)
    (workdir / "status/SUCCESS").write_text("previous run\n", encoding="utf-8")
    (workdir / "status/FAILED").write_text("stale failure\n", encoding="utf-8")

    args = parse_args(
        [
            "run",
            str(workdir),
            "--reference",
            str(reference),
            "--stranded",
            "reverse",
            "--threads",
            "2",
            "--adapter-r1",
            "GGGG",
            "--adapter-r2",
            "TTTT",
            "--dry-run",
        ]
    )
    assert execute(args) == 0

    assert (workdir / "status/manifest.json").read_bytes() == manifest_before
    assert all(step_done_file(workdir, phase).exists() for phase in PHASES if phase != "align")
    assert (workdir / "logs/run_parameters.txt").read_bytes() == authoritative_params
    assert (workdir / "metadata/versions.txt").read_bytes() == authoritative_versions
    assert (workdir / "logs/dry_run_parameters.txt").exists()
    assert (workdir / "metadata/dry_run_versions.txt").exists()
    assert (workdir / "status/DRY_RUN_OK").exists()
    assert not (workdir / "status/SUCCESS").exists()
    assert not (workdir / "status/FAILED").exists()
    assert "would invalidate cutadapt" in capsys.readouterr().out


@pytest.mark.skipif(shutil.which("cutadapt") is None, reason="real Cutadapt is not installed")
def test_real_cutadapt_removes_configured_adapter_and_report_validates(tmp_path):
    input_path = tmp_path / "fastq/sample.fastq.gz"
    input_path.parent.mkdir()
    with gzip.open(input_path, "wt", encoding="utf-8") as f:
        f.write("@read1\nACGTACGTAGATCGGAAGAGC\n+\nIIIIIIIIIIIIIIIIIIIII\n")
    sample = Sample("sample", input_path)
    trimming = TrimmingConfig(TRIMMING_ADAPTERS, ("AGATCGGAAGAGC",), ())
    command = build_cutadapt_command(tmp_path, sample, 2, LAYOUT_SINGLE, trimming)
    (tmp_path / "trimmed").mkdir()
    (tmp_path / "qc/cutadapt").mkdir(parents=True)

    subprocess.run(command, check=True, capture_output=True, text=True)

    validate_cutadapt_report(tmp_path / "qc/cutadapt/sample.cutadapt.json", trimming, LAYOUT_SINGLE)
    with gzip.open(tmp_path / "trimmed/sample.trimmed.fastq.gz", "rt", encoding="utf-8") as f:
        assert f.read().splitlines()[1] == "ACGTACGT"


@pytest.mark.skipif(shutil.which("cutadapt") is None, reason="real Cutadapt is not installed")
def test_real_cutadapt_explicit_no_trimming_preserves_read_and_validates_empty_adapters(tmp_path):
    input_path = tmp_path / "fastq/sample.fastq.gz"
    input_path.parent.mkdir()
    sequence = "ACGTACGTAGATCGGAAGAGC"
    with gzip.open(input_path, "wt", encoding="utf-8") as f:
        f.write(f"@read1\n{sequence}\n+\n{'I' * len(sequence)}\n")
    sample = Sample("sample", input_path)
    trimming = TrimmingConfig(TRIMMING_DISABLED, (), ())
    command = build_cutadapt_command(tmp_path, sample, 2, LAYOUT_SINGLE, trimming)
    (tmp_path / "trimmed").mkdir()
    (tmp_path / "qc/cutadapt").mkdir(parents=True)

    subprocess.run(command, check=True, capture_output=True, text=True)

    validate_cutadapt_report(tmp_path / "qc/cutadapt/sample.cutadapt.json", trimming, LAYOUT_SINGLE)
    with gzip.open(tmp_path / "trimmed/sample.trimmed.fastq.gz", "rt", encoding="utf-8") as f:
        assert f.read().splitlines()[1] == sequence
