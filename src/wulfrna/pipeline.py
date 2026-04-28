from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .io import capture_versions_file, now_iso

PHASES = ["fastqc_raw", "cutadapt", "fastqc_trimmed", "quant", "aggregate", "multiqc"]


class PipelineError(Exception):
    def __init__(self, message: str, step: str, sample: Optional[str] = None):
        super().__init__(message)
        self.step = step
        self.sample = sample


def run_cmd(cmd: List[str], log_path: Path, step: str, sample: Optional[str] = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise PipelineError(f"Command failed ({proc.returncode}): {' '.join(cmd)}", step=step, sample=sample)


def check_tools(tools: List[str]) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise PipelineError(f"Missing required tools in PATH: {', '.join(missing)}", step="tool_check")


def resolve_reference_dir(reference_root: Path, genome: Optional[str]) -> Path:
    if genome is None:
        return reference_root
    reference_dir = reference_root / genome
    if not reference_dir.is_dir():
        raise PipelineError(f"Genome-specific reference directory not found: {reference_dir}", step="reference_check")
    return reference_dir


def validate_reference(reference_dir: Path, quantifier: str) -> Dict[str, Path]:
    tx2gene = reference_dir / "combined_tx2gene.tsv"
    missing: List[str] = []

    refs: Dict[str, Path] = {"tx2gene": tx2gene}
    if quantifier == "salmon":
        salmon_index = reference_dir / "salmon_index"
        refs["index"] = salmon_index
        if not salmon_index.is_dir():
            missing.append(str(salmon_index))
    elif quantifier == "kallisto":
        kallisto_index = reference_dir / "kallisto_index" / "combined_transcripts.kidx"
        refs["index"] = kallisto_index
        if not kallisto_index.is_file():
            missing.append(str(kallisto_index))
    else:
        raise PipelineError(f"Unsupported quantifier: {quantifier}", step="reference_check")

    if not tx2gene.is_file():
        missing.append(str(tx2gene))

    if missing:
        raise PipelineError("Reference directory is missing required files: " + ", ".join(missing), step="reference_check")

    return refs


def detect_samples(workdir: Path) -> List[Tuple[str, Path, Path]]:
    fastq_dir = workdir / "fastq"
    if not fastq_dir.is_dir():
        raise PipelineError("workdir/fastq directory does not exist", step="sample_detection")

    r1_files = sorted(fastq_dir.glob("*_R1.fastq.gz"))
    if not r1_files:
        raise PipelineError("No paired FASTQ files found in fastq/", step="sample_detection")

    samples: Dict[str, Tuple[Path, Path]] = {}
    pattern = re.compile(r"^(?P<sample>.+)_R1\.fastq\.gz$")

    for r1 in r1_files:
        m = pattern.match(r1.name)
        if not m:
            continue
        sample_id = m.group("sample")
        if not sample_id:
            raise PipelineError(f"Invalid sample id from filename: {r1.name}", step="sample_detection")
        if sample_id in samples:
            raise PipelineError(f"Conflicting sample interpretation for sample_id={sample_id}", step="sample_detection", sample=sample_id)
        r2 = fastq_dir / f"{sample_id}_R2.fastq.gz"
        if not r2.exists():
            raise PipelineError(f"Missing R2 FASTQ for sample {sample_id}: {r2}", step="sample_detection", sample=sample_id)
        samples[sample_id] = (r1, r2)

    if not samples:
        raise PipelineError("No valid sample pairs detected", step="sample_detection")

    for r2 in fastq_dir.glob("*_R2.fastq.gz"):
        candidate = r2.name[: -len("_R2.fastq.gz")]
        if candidate not in samples:
            raise PipelineError(f"Orphan R2 FASTQ without matching R1: {r2.name}", step="sample_detection", sample=candidate)

    sample_list = sorted((sample, pair[0], pair[1]) for sample, pair in samples.items())
    samples_tsv = workdir / "samples.tsv"
    with samples_tsv.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample_id", "r1", "r2"])
        for sample_id, r1, r2 in sample_list:
            w.writerow([sample_id, str(r1.relative_to(workdir)), str(r2.relative_to(workdir))])

    return sample_list


def write_params(workdir: Path, args: argparse.Namespace) -> None:
    out = workdir / "logs" / "run_parameters.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        for k, v in vars(args).items():
            f.write(f"{k}: {v}\n")


def capture_versions(workdir: Path) -> None:
    capture_versions_file(workdir / "metadata" / "versions.txt")


def parse_tx2gene(tx2gene_tsv: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with tx2gene_tsv.open("r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        first_row = next(reader, None)
        if first_row is None:
            raise PipelineError(f"Empty tx2gene file: {tx2gene_tsv}", step="tx2gene")

        tx_idx = 0
        gene_idx = 1
        rows = []
        lower_first_row = [h.strip().lower() for h in first_row]
        if "transcript_id" in lower_first_row and "gene_id" in lower_first_row:
            header = first_row
            lower_header = [h.strip().lower() for h in header]
            tx_idx = lower_header.index("transcript_id")
            gene_idx = lower_header.index("gene_id")
        else:
            rows.append((1, first_row))

        for row_num, row in enumerate(reader, start=2):
            rows.append((row_num, row))

        for row_num, row in rows:
            if len(row) <= max(tx_idx, gene_idx):
                raise PipelineError(f"Malformed tx2gene row {row_num} in {tx2gene_tsv}", step="tx2gene")
            tx_id = row[tx_idx].strip()
            gene_id = row[gene_idx].strip()
            if not tx_id or not gene_id:
                raise PipelineError(f"Blank transcript or gene ID at row {row_num} in {tx2gene_tsv}", step="tx2gene")
            mapping[tx_id] = gene_id

    if not mapping:
        raise PipelineError(f"No mappings found in {tx2gene_tsv}", step="tx2gene")
    return mapping


def run_salmon_quant(workdir: Path, sample: str, refs: Dict[str, Path], threads: int, stranded: str) -> Path:
    trimmed_r1 = workdir / "trimmed" / f"{sample}_R1.trimmed.fastq.gz"
    trimmed_r2 = workdir / "trimmed" / f"{sample}_R2.trimmed.fastq.gz"
    sample_dir = workdir / "abundance" / "salmon" / sample
    salmon_libtype = {"none": "IU", "forward": "ISF", "reverse": "ISR"}[stranded]
    run_cmd(
        [
            "salmon",
            "quant",
            "-i",
            str(refs["index"]),
            "-l",
            salmon_libtype,
            "--validateMappings",
            "--seqBias",
            "--gcBias",
            "-1",
            str(trimmed_r1),
            "-2",
            str(trimmed_r2),
            "-p",
            str(threads),
            "-o",
            str(sample_dir),
        ],
        workdir / "logs/steps" / f"{sample}.salmon.log",
        step="salmon_quant",
        sample=sample,
    )
    quant_sf = sample_dir / "quant.sf"
    if not quant_sf.exists():
        raise PipelineError(f"Missing Salmon quant file for sample {sample}: {quant_sf}", step="salmon_quant", sample=sample)
    return quant_sf


def run_kallisto_quant(workdir: Path, sample: str, refs: Dict[str, Path], threads: int, stranded: str) -> Path:
    trimmed_r1 = workdir / "trimmed" / f"{sample}_R1.trimmed.fastq.gz"
    trimmed_r2 = workdir / "trimmed" / f"{sample}_R2.trimmed.fastq.gz"
    sample_dir = workdir / "abundance" / "kallisto" / sample
    cmd = [
        "kallisto",
        "quant",
        "-i",
        str(refs["index"]),
        "-t",
        str(threads),
        "-o",
        str(sample_dir),
    ]
    if stranded == "forward":
        cmd.append("--fr-stranded")
    elif stranded == "reverse":
        cmd.append("--rf-stranded")
    cmd.extend([str(trimmed_r1), str(trimmed_r2)])
    run_cmd(
        cmd,
        workdir / "logs/steps" / f"{sample}.kallisto.log",
        step="kallisto_quant",
        sample=sample,
    )
    abundance_tsv = sample_dir / "abundance.tsv"
    if not abundance_tsv.exists():
        raise PipelineError(
            f"Missing kallisto abundance file for sample {sample}: {abundance_tsv}",
            step="kallisto_quant",
            sample=sample,
        )
    return abundance_tsv


def aggregate_transcript_quant(
    sample_quant_files: Dict[str, Path],
    tx2gene_map: Dict[str, str],
    quantifier: str,
    out_expected: Path,
    out_tpm: Path,
    mapping_stats_out: Path,
    min_mapping_rate: float,
) -> None:
    sample_order = sorted(sample_quant_files)
    expected_by_gene: Dict[str, Dict[str, float]] = defaultdict(dict)
    tpm_by_gene: Dict[str, Dict[str, float]] = defaultdict(dict)

    if quantifier == "salmon":
        tx_col, exp_col, tpm_col = "Name", "NumReads", "TPM"
        step = "aggregate_salmon"
    elif quantifier == "kallisto":
        tx_col, exp_col, tpm_col = "target_id", "est_counts", "tpm"
        step = "aggregate_kallisto"
    else:
        raise PipelineError(f"Unsupported quantifier for aggregation: {quantifier}", step="aggregate")

    mapping_stats_out.parent.mkdir(parents=True, exist_ok=True)
    with mapping_stats_out.open("w", encoding="utf-8", newline="") as stats_f:
        stats_writer = csv.writer(stats_f, delimiter="\t")
        stats_writer.writerow(["sample_id", "total_transcripts", "mapped_transcripts", "mapping_rate"])
        for sample in sample_order:
            sample_exp: Dict[str, float] = defaultdict(float)
            sample_tpm: Dict[str, float] = defaultdict(float)

            path = sample_quant_files[sample]
            total_transcripts = 0
            mapped_transcripts = 0
            with path.open("r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                if reader.fieldnames is None:
                    raise PipelineError(f"Quantification file is missing header: {path}", step=step, sample=sample)

                required_cols = {tx_col, exp_col, tpm_col}
                if not required_cols.issubset(set(reader.fieldnames)):
                    raise PipelineError(
                        f"Quantification file missing required columns ({', '.join(sorted(required_cols))}): {path}",
                        step=step,
                        sample=sample,
                    )

                for row in reader:
                    total_transcripts += 1
                    tx_id = row[tx_col]
                    gene_id = tx2gene_map.get(tx_id)
                    if gene_id is None:
                        continue
                    mapped_transcripts += 1
                    sample_exp[gene_id] += float(row[exp_col])
                    sample_tpm[gene_id] += float(row[tpm_col])

            if total_transcripts == 0:
                raise PipelineError(f"No transcript records found in quantification file: {path}", step=step, sample=sample)

            mapping_rate = mapped_transcripts / total_transcripts
            stats_writer.writerow([sample, total_transcripts, mapped_transcripts, f"{mapping_rate:.6f}"])
            if mapping_rate < min_mapping_rate:
                raise PipelineError(
                    (
                        f"Low tx2gene mapping rate for sample {sample}: "
                        f"{mapping_rate:.4f} < {min_mapping_rate:.4f} "
                        f"(mapped {mapped_transcripts}/{total_transcripts})"
                    ),
                    step=step,
                    sample=sample,
                )

            for gene_id, value in sample_exp.items():
                expected_by_gene[gene_id][sample] = value
            for gene_id, value in sample_tpm.items():
                tpm_by_gene[gene_id][sample] = value

    write_gene_matrix(out_expected, sample_order, expected_by_gene)
    write_gene_matrix(out_tpm, sample_order, tpm_by_gene)


def write_gene_matrix(out_file: Path, sample_order: List[str], matrix: Dict[str, Dict[str, float]]) -> None:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    with out_file.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene_id", *sample_order])
        for gid in sorted(matrix):
            values = [f"{matrix[gid].get(sample, 0.0):.6f}" for sample in sample_order]
            w.writerow([gid, *values])


def write_summary(
    workdir: Path,
    status: str,
    started: str,
    finished: str,
    sample_count: int,
    completed_count: int,
    failed_count: int,
) -> None:
    status_dir = workdir / "status"
    out = status_dir / "summary.txt"
    with out.open("w", encoding="utf-8") as f:
        f.write(f"Pipeline status: {status}\n")
        f.write(f"Started: {started}\n")
        f.write(f"Finished: {finished}\n")
        f.write(f"Samples detected: {sample_count}\n")
        f.write(f"Samples completed: {completed_count}\n")
        f.write(f"Failed samples: {failed_count}\n")
        if status == "DRY_RUN_OK":
            f.write("Expected outputs (not generated in dry-run):\n")
        else:
            f.write("Main outputs:\n")
        f.write("- abundance/gene_expected_counts.tsv\n")
        f.write("- abundance/gene_tpm.tsv\n")
        f.write("- multiqc/multiqc_report.html\n")


def ensure_dirs(workdir: Path) -> None:
    required = [
        "trimmed",
        "qc/fastqc_raw",
        "qc/fastqc_trimmed",
        "qc/cutadapt",
        "abundance",
        "abundance/salmon",
        "abundance/kallisto",
        "multiqc",
        "logs",
        "logs/steps",
        "status",
        "status/steps",
        "metadata",
    ]
    for d in required:
        (workdir / d).mkdir(parents=True, exist_ok=True)


def phase_index(phase: str) -> int:
    if phase not in PHASES:
        raise PipelineError(f"Unknown phase: {phase}", step="resume")
    return PHASES.index(phase)


def step_done_file(workdir: Path, phase: str) -> Path:
    return workdir / "status" / "steps" / f"{phase}.done"


def mark_phase_done(workdir: Path, phase: str) -> None:
    step_done_file(workdir, phase).write_text(now_iso() + "\n", encoding="utf-8")


def output_exists_nonempty(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def fastqc_prefix(fastq_path: Path) -> str:
    name = fastq_path.name
    for suffix in [".fastq.gz", ".fq.gz", ".fastq", ".fq"]:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return fastq_path.stem


def phase_outputs(workdir: Path, phase: str, samples: List[Tuple[str, Path, Path]], quantifier: str) -> List[Path]:
    outputs: List[Path] = []
    if phase == "fastqc_raw":
        outdir = workdir / "qc" / "fastqc_raw"
        for _, r1, r2 in samples:
            for fq in [r1, r2]:
                prefix = fastqc_prefix(fq)
                outputs.append(outdir / f"{prefix}_fastqc.html")
                outputs.append(outdir / f"{prefix}_fastqc.zip")
    elif phase == "cutadapt":
        for sample, _, _ in samples:
            outputs.extend(
                [
                    workdir / "trimmed" / f"{sample}_R1.trimmed.fastq.gz",
                    workdir / "trimmed" / f"{sample}_R2.trimmed.fastq.gz",
                    workdir / "qc" / "cutadapt" / f"{sample}.cutadapt.json",
                ]
            )
    elif phase == "fastqc_trimmed":
        outdir = workdir / "qc" / "fastqc_trimmed"
        for sample, _, _ in samples:
            for fq in [workdir / "trimmed" / f"{sample}_R1.trimmed.fastq.gz", workdir / "trimmed" / f"{sample}_R2.trimmed.fastq.gz"]:
                prefix = fastqc_prefix(fq)
                outputs.append(outdir / f"{prefix}_fastqc.html")
                outputs.append(outdir / f"{prefix}_fastqc.zip")
    elif phase == "quant":
        if quantifier == "salmon":
            for sample, _, _ in samples:
                outputs.append(workdir / "abundance" / "salmon" / sample / "quant.sf")
        else:
            for sample, _, _ in samples:
                outputs.append(workdir / "abundance" / "kallisto" / sample / "abundance.tsv")
    elif phase == "aggregate":
        outputs.extend([workdir / "abundance/gene_expected_counts.tsv", workdir / "abundance/gene_tpm.tsv"])
    elif phase == "multiqc":
        outputs.append(workdir / "multiqc/multiqc_report.html")
    else:
        raise PipelineError(f"Unknown phase: {phase}", step="resume")
    return outputs


def outputs_exist(paths: List[Path]) -> bool:
    return all(output_exists_nonempty(path) for path in paths)


def is_phase_complete(workdir: Path, phase: str, samples: List[Tuple[str, Path, Path]], quantifier: str) -> bool:
    return step_done_file(workdir, phase).exists() and outputs_exist(phase_outputs(workdir, phase, samples, quantifier))


def fingerprint_file(path: Path) -> Dict[str, str]:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    stat = path.stat()
    return {"path": str(path), "size": str(stat.st_size), "mtime_ns": str(stat.st_mtime_ns), "sha256": h.hexdigest()}


def load_manifest(workdir: Path) -> Optional[Dict[str, object]]:
    manifest_path = workdir / "status" / "manifest.json"
    if not manifest_path.exists():
        return None
    with manifest_path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_manifest(workdir: Path, manifest: Dict[str, object]) -> None:
    manifest_path = workdir / "status" / "manifest.json"
    existing = load_manifest(workdir)
    if existing is not None and "created_at" in existing:
        manifest["created_at"] = existing["created_at"]
    elif "created_at" not in manifest:
        manifest["created_at"] = now_iso()
    manifest["updated_at"] = now_iso()
    with manifest_path.open("w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")


def compare_manifest(existing: Dict[str, object], current: Dict[str, object]) -> Tuple[Optional[str], List[str]]:
    existing_samples = sorted(existing.get("sample_ids", []))
    current_samples = sorted(current.get("sample_ids", []))
    if existing_samples != current_samples:
        raise PipelineError(
            "Sample set changed since last run; refusing automatic resume. Adjust inputs or run in a new workdir.",
            step="resume",
        )

    forced_phase: Optional[str] = None
    reasons: List[str] = []

    def force_from(phase: str, reason: str) -> None:
        nonlocal forced_phase
        reasons.append(reason)
        if forced_phase is None or phase_index(phase) < phase_index(forced_phase):
            forced_phase = phase

    if existing.get("quantifier") != current.get("quantifier"):
        force_from("quant", "quantifier changed")
    if existing.get("reference_dir") != current.get("reference_dir"):
        force_from("quant", "reference_dir changed")
    if existing.get("stranded") != current.get("stranded"):
        force_from("quant", "stranded changed")

    old_fp = existing.get("tx2gene_fingerprint", {})
    new_fp = current.get("tx2gene_fingerprint", {})
    if old_fp != new_fp:
        force_from("aggregate", "combined_tx2gene.tsv fingerprint changed")

    return forced_phase, reasons


def should_run_phase(
    workdir: Path,
    phase: str,
    samples: List[Tuple[str, Path, Path]],
    quantifier: str,
    no_resume: bool,
    effective_force_from: Optional[str],
) -> bool:
    if no_resume:
        return True
    if effective_force_from is not None and phase_index(phase) >= phase_index(effective_force_from):
        return True
    return not is_phase_complete(workdir, phase, samples, quantifier)


def execute(args: argparse.Namespace) -> int:
    workdir = Path(args.workdir).resolve()
    reference_root = Path(args.reference).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    ensure_dirs(workdir)

    status_dir = workdir / "status"
    running = status_dir / "RUNNING"
    running.write_text(now_iso() + "\n", encoding="utf-8")
    started = now_iso()

    completed_samples = 0
    samples: List[Tuple[str, Path, Path]] = []

    try:
        base_tools = ["fastqc", "cutadapt", "multiqc"]
        quant_tool = "salmon" if args.quantifier == "salmon" else "kallisto"
        check_tools(base_tools + [quant_tool])

        reference_dir = resolve_reference_dir(reference_root, args.genome)
        refs = validate_reference(reference_dir, args.quantifier)
        samples = detect_samples(workdir)
        tx2gene_map = parse_tx2gene(refs["tx2gene"])
        write_params(workdir, args)
        capture_versions(workdir)

        manifest: Dict[str, object] = {
            "workdir": str(workdir),
            "reference_dir": str(reference_dir),
            "quantifier": args.quantifier,
            "stranded": args.stranded,
            "sample_ids": [sample for sample, _, _ in samples],
            "reference_files": {
                "index": str(refs["index"]),
                "combined_tx2gene_tsv": str(refs["tx2gene"]),
            },
            "tx2gene_fingerprint": fingerprint_file(refs["tx2gene"]),
        }
        existing_manifest = load_manifest(workdir)
        auto_force_from: Optional[str] = None
        manifest_reasons: List[str] = []
        if existing_manifest is not None:
            auto_force_from, manifest_reasons = compare_manifest(existing_manifest, manifest)
            if auto_force_from is not None:
                print(f"[force] {auto_force_from} (manifest: {', '.join(manifest_reasons)})")
        write_manifest(workdir, manifest)

        if args.dry_run:
            (status_dir / "DRY_RUN_OK").write_text(now_iso() + "\n", encoding="utf-8")
            (status_dir / "finished_at.txt").write_text(now_iso() + "\n", encoding="utf-8")
            if running.exists():
                running.unlink()
            write_summary(
                workdir,
                "DRY_RUN_OK",
                started,
                now_iso(),
                len(samples),
                0,
                0,
            )
            return 0

        effective_force_from = args.force_from
        if auto_force_from is not None and (
            effective_force_from is None or phase_index(auto_force_from) < phase_index(effective_force_from)
        ):
            effective_force_from = auto_force_from

        quant_map: Dict[str, Path] = {}
        for phase in PHASES:
            phase_should_run = should_run_phase(workdir, phase, samples, args.quantifier, args.no_resume, effective_force_from)
            if not phase_should_run:
                print(f"[skip] {phase}")
                continue
            if args.no_resume:
                print(f"[resume disabled] {phase}")
            elif effective_force_from is not None and phase_index(phase) >= phase_index(effective_force_from):
                print(f"[force] {phase}")
            else:
                print(f"[run] {phase}")

            if phase == "fastqc_raw":
                for sample, r1, r2 in samples:
                    run_cmd(
                        ["fastqc", "-t", str(max(1, args.threads // 2)), "-o", str(workdir / "qc/fastqc_raw"), str(r1), str(r2)],
                        workdir / "logs/steps" / f"{sample}.fastqc_raw.log",
                        step="fastqc_raw",
                        sample=sample,
                    )
            elif phase == "cutadapt":
                for sample, r1, r2 in samples:
                    trimmed_r1 = workdir / "trimmed" / f"{sample}_R1.trimmed.fastq.gz"
                    trimmed_r2 = workdir / "trimmed" / f"{sample}_R2.trimmed.fastq.gz"
                    run_cmd(
                        [
                            "cutadapt",
                            "-j",
                            str(max(1, args.threads // 2)),
                            "-o",
                            str(trimmed_r1),
                            "-p",
                            str(trimmed_r2),
                            "--json",
                            str(workdir / "qc/cutadapt" / f"{sample}.cutadapt.json"),
                            str(r1),
                            str(r2),
                        ],
                        workdir / "logs/steps" / f"{sample}.cutadapt.log",
                        step="cutadapt",
                        sample=sample,
                    )
            elif phase == "fastqc_trimmed":
                for sample, _, _ in samples:
                    trimmed_r1 = workdir / "trimmed" / f"{sample}_R1.trimmed.fastq.gz"
                    trimmed_r2 = workdir / "trimmed" / f"{sample}_R2.trimmed.fastq.gz"
                    run_cmd(
                        ["fastqc", "-t", str(max(1, args.threads // 2)), "-o", str(workdir / "qc/fastqc_trimmed"), str(trimmed_r1), str(trimmed_r2)],
                        workdir / "logs/steps" / f"{sample}.fastqc_trimmed.log",
                        step="fastqc_trimmed",
                        sample=sample,
                    )
            elif phase == "quant":
                for sample, _, _ in samples:
                    if args.quantifier == "salmon":
                        quant_map[sample] = run_salmon_quant(workdir, sample, refs, args.threads, args.stranded)
                    else:
                        quant_map[sample] = run_kallisto_quant(workdir, sample, refs, args.threads, args.stranded)
                completed_samples = len(samples)
            elif phase == "aggregate":
                if not quant_map:
                    for sample, _, _ in samples:
                        if args.quantifier == "salmon":
                            quant_map[sample] = workdir / "abundance" / "salmon" / sample / "quant.sf"
                        else:
                            quant_map[sample] = workdir / "abundance" / "kallisto" / sample / "abundance.tsv"
                aggregate_transcript_quant(
                    quant_map,
                    tx2gene_map,
                    args.quantifier,
                    workdir / "abundance/gene_expected_counts.tsv",
                    workdir / "abundance/gene_tpm.tsv",
                    workdir / "logs" / "tx2gene_mapping_stats.tsv",
                    args.min_mapping_rate,
                )
            elif phase == "multiqc":
                run_cmd(
                    ["multiqc", str(workdir), "-o", str(workdir / "multiqc")],
                    workdir / "logs/steps" / "multiqc.log",
                    step="multiqc",
                    sample=None,
                )
            else:
                raise PipelineError(f"Unhandled phase: {phase}", step="resume")

            if not outputs_exist(phase_outputs(workdir, phase, samples, args.quantifier)):
                raise PipelineError(f"Expected outputs missing after phase: {phase}", step=phase)
            mark_phase_done(workdir, phase)
            print(f"[done] {phase}")

        required_outputs = [
            workdir / "abundance/gene_expected_counts.tsv",
            workdir / "abundance/gene_tpm.tsv",
            workdir / "multiqc/multiqc_report.html",
        ]
        for path in required_outputs:
            if not path.exists():
                raise PipelineError(f"Required output missing: {path}", step="final_validation")

        (status_dir / "SUCCESS").write_text(now_iso() + "\n", encoding="utf-8")
        (status_dir / "finished_at.txt").write_text(now_iso() + "\n", encoding="utf-8")
        if running.exists():
            running.unlink()
        write_manifest(workdir, manifest)

        write_summary(
            workdir,
            "SUCCESS",
            started,
            now_iso(),
            len(samples),
            completed_samples,
            0,
        )
        return 0

    except PipelineError as e:
        (status_dir / "FAILED").write_text(now_iso() + "\n", encoding="utf-8")
        (status_dir / "failed_step.txt").write_text(e.step + "\n", encoding="utf-8")
        if e.sample is not None:
            (status_dir / "failed_sample.txt").write_text(e.sample + "\n", encoding="utf-8")

        if running.exists():
            running.unlink()

        write_summary(
            workdir,
            "FAILED",
            started,
            now_iso(),
            len(samples),
            completed_samples,
            max(0, len(samples) - completed_samples),
        )
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
