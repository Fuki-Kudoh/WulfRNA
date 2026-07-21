import stat
from pathlib import Path

import pytest

from wulfrna.cli import parse_args
from wulfrna.pipeline import STAR_INDEX_REQUIRED_FILES, execute


def write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def install_fake_tools(bin_dir: Path, include_star: bool = True) -> None:
    bin_dir.mkdir()
    write_executable(
        bin_dir / "fastqc",
        r"""#!/usr/bin/env python3
from pathlib import Path
import sys
args = sys.argv[1:]
outdir = Path(args[args.index('-o') + 1])
outdir.mkdir(parents=True, exist_ok=True)
for item in args:
    if item.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')):
        name = Path(item).name
        for suffix in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
            if name.endswith(suffix):
                prefix = name[:-len(suffix)]
                break
        (outdir / f'{prefix}_fastqc.html').write_text('html\n', encoding='utf-8')
        (outdir / f'{prefix}_fastqc.zip').write_text('zip\n', encoding='utf-8')
""",
    )
    write_executable(
        bin_dir / "cutadapt",
        r"""#!/usr/bin/env python3
from pathlib import Path
import sys
args = sys.argv[1:]
outputs = []
for flag in ['-o', '-p']:
    if flag in args:
        outputs.append(Path(args[args.index(flag) + 1]))
if '--json' in args:
    outputs.append(Path(args[args.index('--json') + 1]))
for output in outputs:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text('cutadapt\n', encoding='utf-8')
""",
    )
    write_executable(
        bin_dir / "salmon",
        r"""#!/usr/bin/env python3
from pathlib import Path
import sys
args = sys.argv[1:]
outdir = Path(args[args.index('-o') + 1])
outdir.mkdir(parents=True, exist_ok=True)
(outdir / 'quant.sf').write_text('Name\tLength\tEffectiveLength\tTPM\tNumReads\ntx1\t100\t100\t1.0\t5\n', encoding='utf-8')
""",
    )
    write_executable(
        bin_dir / "kallisto",
        r"""#!/usr/bin/env python3
from pathlib import Path
import sys
args = sys.argv[1:]
outdir = Path(args[args.index('-o') + 1])
outdir.mkdir(parents=True, exist_ok=True)
(outdir / 'abundance.tsv').write_text('target_id\tlength\teff_length\test_counts\ttpm\ntx1\t100\t100\t5\t1.0\n', encoding='utf-8')
""",
    )
    write_executable(
        bin_dir / "multiqc",
        r"""#!/usr/bin/env python3
from pathlib import Path
import sys
args = sys.argv[1:]
outdir = Path(args[args.index('-o') + 1])
outdir.mkdir(parents=True, exist_ok=True)
(outdir / 'multiqc_report.html').write_text('multiqc\n', encoding='utf-8')
""",
    )
    if include_star:
        write_executable(
            bin_dir / "STAR",
            r"""#!/usr/bin/env python3
from pathlib import Path
import sys
args = sys.argv[1:]
if '--version' in args:
    print('STAR 2.7.fake')
    raise SystemExit(0)
prefix = Path(args[args.index('--outFileNamePrefix') + 1])
prefix.mkdir(parents=True, exist_ok=True)
for name in ['Aligned.sortedByCoord.out.bam', 'SJ.out.tab', 'Log.final.out']:
    (prefix / name).write_text(f'{name}\n', encoding='utf-8')
(prefix / 'ReadsPerGene.out.tab').write_text('N_unmapped\t0\t0\t0\nN_multimapping\t0\t0\t0\nN_noFeature\t0\t0\t0\nN_ambiguous\t0\t0\t0\ngene1\t5\t6\t7\n', encoding='utf-8')
""",
        )
        write_executable(
            bin_dir / "samtools",
            r"""#!/usr/bin/env python3
from pathlib import Path
import sys
if '--version' in sys.argv[1:]:
    print('samtools 1.fake')
    raise SystemExit(0)
bam = Path(sys.argv[2])
Path(str(bam) + '.bai').write_text('bai\n', encoding='utf-8')
""",
        )


def make_reference(ref_dir: Path) -> None:
    (ref_dir / "combined_tx2gene.tsv").write_text("tx1\tgene1\n", encoding="utf-8")
    (ref_dir / "combined_gene_annotation.tsv").write_text("gene_id\tGeneName\tgene_length_bp\ngene1\tGene1\t1000\n", encoding="utf-8")
    (ref_dir / "salmon_index").mkdir()
    (ref_dir / "kallisto_index").mkdir()
    (ref_dir / "kallisto_index" / "combined_transcripts.kidx").write_text("index\n", encoding="utf-8")
    star_index = ref_dir / "star_index"
    star_index.mkdir()
    for name in STAR_INDEX_REQUIRED_FILES:
        (star_index / name).write_text(f"{name}\n", encoding="utf-8")


def run_pipeline(tmp_path, monkeypatch, single_end: bool, quantifier: str = "salmon"):
    bin_dir = tmp_path / "bin"
    workdir = tmp_path / "work"
    ref_dir = tmp_path / "ref"
    fastq_dir = workdir / "fastq"
    ref_dir.mkdir()
    fastq_dir.mkdir(parents=True)
    install_fake_tools(bin_dir)
    make_reference(ref_dir)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path('/usr/bin')}")

    if single_end:
        (fastq_dir / "sample.fastq.gz").write_text("reads\n", encoding="utf-8")
        layout_args = ["--single-end"]
    else:
        (fastq_dir / "sample_R1.fastq.gz").write_text("reads\n", encoding="utf-8")
        (fastq_dir / "sample_R2.fastq.gz").write_text("reads\n", encoding="utf-8")
        layout_args = []

    args = parse_args([
        "run",
        str(workdir),
        "--reference",
        str(ref_dir),
        "--stranded",
        "reverse",
        "--threads",
        "2",
        "--aligner",
        "star",
        "--quantifier",
        quantifier,
        *layout_args,
    ])

    assert execute(args) == 0
    return workdir


def assert_star_outputs(workdir: Path) -> None:
    star_dir = workdir / "align" / "star" / "sample"
    for name in [
        "Aligned.sortedByCoord.out.bam",
        "Aligned.sortedByCoord.out.bam.bai",
        "SJ.out.tab",
        "ReadsPerGene.out.tab",
        "Log.final.out",
    ]:
        assert (star_dir / name).stat().st_size > 0
    assert (workdir / "abundance" / "gene_expected_counts.tsv").stat().st_size > 0
    assert (workdir / "abundance" / "gene_tpm.tsv").stat().st_size > 0
    assert (workdir / "abundance" / "salmon_gene_expected_counts.tsv").read_text().splitlines()[0] == "gene_id\tGeneName\tsample"
    assert (workdir / "abundance" / "salmon_gene_tpm.tsv").read_text().splitlines()[0] == "gene_id\tGeneName\tsample"
    assert (workdir / "abundance" / "gene_expected_counts.tsv").read_text().splitlines()[0] == "gene_id\tsample"
    copied = workdir / "abundance" / "star_gene_counts" / "sample.star.ReadsPerGene.out.tab"
    original = workdir / "align" / "star" / "sample" / "ReadsPerGene.out.tab"
    assert copied.read_bytes() == original.read_bytes()
    counts = (workdir / "abundance" / "star_gene_counts.tsv").read_text().splitlines()
    assert counts == ["gene_id\tGeneName\tsample", "gene1\tGene1\t7"]
    tpm = (workdir / "abundance" / "star_gene_tpm.tsv").read_text().splitlines()
    assert tpm == ["gene_id\tGeneName\tsample", "gene1\tGene1\t1000000.000000"]

    star_log = (workdir / "logs" / "steps" / "sample.star.log").read_text(encoding="utf-8")
    assert "--twopassMode Basic" in star_log
    assert "--outSAMattributes NH HI AS nM MD" in star_log
    assert "--outSAMtype BAM SortedByCoordinate" in star_log
    versions = (workdir / "metadata" / "versions.txt").read_text(encoding="utf-8")
    assert "STAR: STAR 2.7.fake" in versions
    assert "samtools: samtools 1.fake" in versions


def test_star_alignment_outputs_for_paired_end(monkeypatch, tmp_path):
    workdir = run_pipeline(tmp_path, monkeypatch, single_end=False)

    assert_star_outputs(workdir)


def test_star_alignment_outputs_for_single_end(monkeypatch, tmp_path):
    workdir = run_pipeline(tmp_path, monkeypatch, single_end=True)

    assert_star_outputs(workdir)


def test_star_tools_not_required_when_aligner_omitted(monkeypatch, tmp_path):
    bin_dir = tmp_path / "bin"
    workdir = tmp_path / "work"
    ref_dir = tmp_path / "ref"
    fastq_dir = workdir / "fastq"
    ref_dir.mkdir()
    fastq_dir.mkdir(parents=True)
    install_fake_tools(bin_dir, include_star=False)
    make_reference(ref_dir)
    (fastq_dir / "sample_R1.fastq.gz").write_text("reads\n", encoding="utf-8")
    (fastq_dir / "sample_R2.fastq.gz").write_text("reads\n", encoding="utf-8")
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path('/usr/bin')}")

    args = parse_args([
        "run",
        str(workdir),
        "--reference",
        str(ref_dir),
        "--stranded",
        "reverse",
        "--threads",
        "2",
        "--dry-run",
    ])

    assert execute(args) == 0


def test_star_tools_required_when_aligner_star(monkeypatch, tmp_path):
    bin_dir = tmp_path / "bin"
    workdir = tmp_path / "work"
    ref_dir = tmp_path / "ref"
    fastq_dir = workdir / "fastq"
    ref_dir.mkdir()
    fastq_dir.mkdir(parents=True)
    install_fake_tools(bin_dir, include_star=False)
    make_reference(ref_dir)
    (fastq_dir / "sample_R1.fastq.gz").write_text("reads\n", encoding="utf-8")
    (fastq_dir / "sample_R2.fastq.gz").write_text("reads\n", encoding="utf-8")
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path('/usr/bin')}")

    args = parse_args([
        "run",
        str(workdir),
        "--reference",
        str(ref_dir),
        "--stranded",
        "reverse",
        "--threads",
        "2",
        "--aligner",
        "star",
        "--dry-run",
    ])

    assert execute(args) == 1
    assert (workdir / "status" / "failed_step.txt").read_text(encoding="utf-8") == "tool_check\n"


@pytest.mark.parametrize("quantifier", ["salmon", "kallisto"])
def test_v021_migration_and_annotation_resume_behavior(monkeypatch, tmp_path, quantifier):
    import json
    import wulfrna.pipeline as pipeline

    workdir = run_pipeline(tmp_path, monkeypatch, single_end=False, quantifier=quantifier)
    ref_dir = tmp_path / "ref"
    args = parse_args([
        "run", str(workdir), "--reference", str(ref_dir),
        "--stranded", "reverse", "--threads", "2", "--aligner", "star", "--quantifier", quantifier,
    ])

    # Model a completed v0.2.1 directory: phase markers and legacy matrices
    # exist, but source-qualified/STAR aggregate products and annotation
    # manifest fields do not.
    for path in [
        workdir / "abundance" / f"{quantifier}_gene_expected_counts.tsv",
        workdir / "abundance" / f"{quantifier}_gene_tpm.tsv",
        workdir / "abundance/star_gene_counts.tsv",
        workdir / "abundance/star_gene_tpm.tsv",
        workdir / "abundance/star_gene_counts/sample.star.ReadsPerGene.out.tab",
    ]:
        path.unlink()
    manifest_path = workdir / "status/manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.pop("gene_annotation_fingerprint")
    manifest["reference_files"].pop("combined_gene_annotation_tsv")
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    def unexpected(*args, **kwargs):
        raise AssertionError("alignment or transcript quantification was unexpectedly rerun")

    monkeypatch.setattr(pipeline, "run_star_align", unexpected)
    monkeypatch.setattr(pipeline, "run_salmon_quant", unexpected)
    monkeypatch.setattr(pipeline, "run_kallisto_quant", unexpected)
    assert execute(args) == 0
    assert (workdir / "abundance/star_gene_counts/sample.star.ReadsPerGene.out.tab").is_file()

    # Annotation content changes invalidate aggregate and downstream only.
    annotation = ref_dir / "combined_gene_annotation.tsv"
    annotation.write_text("gene_id\tGeneName\tgene_length_bp\ngene1\tChanged\t2000\n", encoding="utf-8")
    assert execute(args) == 0
    assert "Changed" in (workdir / "abundance/star_gene_tpm.tsv").read_text(encoding="utf-8")

    # With all outputs and fingerprints unchanged, aggregate is skipped.
    monkeypatch.setattr(pipeline, "aggregate_transcript_quant", unexpected)
    monkeypatch.setattr(pipeline, "aggregate_star_counts", unexpected)
    assert execute(args) == 0


def test_star_copied_files_are_required_for_aggregate_completion(tmp_path):
    from wulfrna.pipeline import LAYOUT_PAIRED, Sample, phase_outputs

    sample = Sample("sample", tmp_path / "sample_R1.fastq.gz", tmp_path / "sample_R2.fastq.gz")
    outputs = phase_outputs(tmp_path, "aggregate", [sample], "salmon", LAYOUT_PAIRED, "star")
    assert tmp_path / "abundance/star_gene_counts/sample.star.ReadsPerGene.out.tab" in outputs
