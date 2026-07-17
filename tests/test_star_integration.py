import stat
from pathlib import Path

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
prefix = Path(args[args.index('--outFileNamePrefix') + 1])
prefix.mkdir(parents=True, exist_ok=True)
for name in ['Aligned.sortedByCoord.out.bam', 'SJ.out.tab', 'ReadsPerGene.out.tab', 'Log.final.out']:
    (prefix / name).write_text(f'{name}\n', encoding='utf-8')
""",
        )
        write_executable(
            bin_dir / "samtools",
            r"""#!/usr/bin/env python3
from pathlib import Path
import sys
bam = Path(sys.argv[2])
Path(str(bam) + '.bai').write_text('bai\n', encoding='utf-8')
""",
        )


def make_reference(ref_dir: Path) -> None:
    (ref_dir / "combined_tx2gene.tsv").write_text("tx1\tgene1\n", encoding="utf-8")
    (ref_dir / "salmon_index").mkdir()
    star_index = ref_dir / "star_index"
    star_index.mkdir()
    for name in STAR_INDEX_REQUIRED_FILES:
        (star_index / name).write_text(f"{name}\n", encoding="utf-8")


def run_pipeline(tmp_path, monkeypatch, single_end: bool):
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
