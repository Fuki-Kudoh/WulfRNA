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
    run_parser.add_argument("--threads", required=True, type=int, help="Total threads")
    run_parser.add_argument("--quantifier", choices=["salmon", "kallisto"], default="salmon", help="Transcript quantification backend")
    run_parser.add_argument("--dry-run", action="store_true", help="Validate inputs, write metadata, and exit without running analysis tools")

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
    return args


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)
    return execute(args)


if __name__ == "__main__":
    sys.exit(main())
