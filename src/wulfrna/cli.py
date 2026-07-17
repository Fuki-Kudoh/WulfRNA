from __future__ import annotations

import argparse
import sys
from typing import List, Optional

from .pipeline import execute


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="WulfRNA pipeline runner")
    subparsers = parser.add_subparsers(dest="command")

    run_parser = subparsers.add_parser("run", help="Run the RNA-seq pipeline")
    run_parser.add_argument("workdir", help="Working directory containing fastq/")
    run_parser.add_argument("--reference", required=True, help="Reference root directory")
    run_parser.add_argument("--genome", required=False, help="Genome key under reference root (optional)")
    run_parser.add_argument("--stranded", required=True, choices=["none", "forward", "reverse"], help="Library strandedness")
    run_parser.add_argument(
        "--min-mapping-rate",
        required=False,
        type=float,
        default=0.90,
        help="Minimum acceptable tx2gene transcript mapping rate per sample (0.0-1.0)",
    )
    run_parser.add_argument("--threads", required=True, type=int, help="Total threads")
    run_parser.add_argument("--quantifier", choices=["salmon", "kallisto"], default="salmon", help="Transcript quantification backend")
    run_parser.add_argument(
        "--aligner",
        choices=["none", "star"],
        default="none",
        help="Optional alignment backend; 'star' runs STAR alignment and writes isolated outputs under align/star/<sample>/",
    )
    run_parser.add_argument("--single-end", "--SE", dest="single_end", action="store_true", help="Use single-end FASTQ input layout")
    run_parser.add_argument("--fragment-length", type=float, help="Single-end fragment length for kallisto")
    run_parser.add_argument("--fragment-sd", type=float, help="Single-end fragment length standard deviation for kallisto")
    run_parser.add_argument("--dry-run", action="store_true", help="Validate inputs, write metadata, and exit without running analysis tools")
    run_parser.add_argument("--no-resume", action="store_true", help="Disable automatic phase-level resume and rerun all phases")
    run_parser.add_argument(
        "--force-from",
        choices=["fastqc_raw", "cutadapt", "fastqc_trimmed", "align", "quant", "aggregate", "multiqc"],
        help="Force rerun from the selected phase onward",
    )

    return parser


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    args_list = list(sys.argv[1:] if argv is None else argv)
    # Backward compatibility: allow legacy invocation without subcommand.
    if args_list and args_list[0] != "run":
        args_list = ["run", *args_list]

    parser = build_parser()
    args = parser.parse_args(args_list)

    if args.command != "run":
        parser.print_help()
        parser.exit(2)
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if args.min_mapping_rate < 0.0 or args.min_mapping_rate > 1.0:
        parser.error("--min-mapping-rate must be between 0.0 and 1.0")
    if args.single_end and args.quantifier == "kallisto":
        if args.fragment_length is None or args.fragment_sd is None:
            parser.error("--single-end --quantifier kallisto requires --fragment-length and --fragment-sd")
        if args.fragment_length <= 0 or args.fragment_sd <= 0:
            parser.error("--fragment-length and --fragment-sd must be > 0")
    return args


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    return execute(args)


if __name__ == "__main__":
    sys.exit(main())
