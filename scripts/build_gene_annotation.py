#!/usr/bin/env python3
"""Build WulfRNA's gene annotation table from exon records in a GTF file."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


ATTRIBUTE_RE = re.compile(r'^\s*([^\s]+)\s+"([^"]*)"\s*$')


class AnnotationBuildError(Exception):
    """Raised when the source GTF cannot produce an unambiguous annotation."""


@dataclass
class GeneRecords:
    chromosome: str
    gene_name: Optional[str]
    intervals: List[Tuple[int, int]]


def parse_attributes(text: str, line_number: int) -> Dict[str, str]:
    attributes: Dict[str, str] = {}
    for item in text.split(";"):
        item = item.strip()
        if not item:
            continue
        match = ATTRIBUTE_RE.match(item)
        if match is None:
            raise AnnotationBuildError(f"Malformed GTF attribute at line {line_number}: {item!r}")
        key, value = match.groups()
        attributes[key] = value
    return attributes


def merge_interval_length(intervals: Iterable[Tuple[int, int]]) -> int:
    """Return the union length of one-based, closed intervals."""
    merged_length = 0
    current_start: Optional[int] = None
    current_end: Optional[int] = None
    for start, end in sorted(intervals):
        if current_start is None:
            current_start, current_end = start, end
        elif start <= current_end + 1:
            current_end = max(current_end, end)
        else:
            merged_length += current_end - current_start + 1
            current_start, current_end = start, end
    if current_start is not None and current_end is not None:
        merged_length += current_end - current_start + 1
    return merged_length


def build_gene_annotation(gtf_path: Path, output_path: Path) -> None:
    genes: Dict[str, GeneRecords] = {}
    with gtf_path.open("r", encoding="utf-8") as gtf:
        for line_number, line in enumerate(gtf, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n\r").split("\t")
            if len(fields) != 9:
                raise AnnotationBuildError(f"GTF line {line_number} must contain exactly nine tab-separated columns")
            chromosome, _, feature, raw_start, raw_end, _, _, _, raw_attributes = fields
            if feature != "exon":
                continue
            try:
                start, end = int(raw_start), int(raw_end)
            except ValueError:
                raise AnnotationBuildError(f"Non-integer exon coordinates at GTF line {line_number}") from None
            if start < 1 or end < start:
                raise AnnotationBuildError(f"Invalid exon interval {raw_start}-{raw_end} at GTF line {line_number}")
            attributes = parse_attributes(raw_attributes, line_number)
            gene_id = attributes.get("gene_id", "")
            if not gene_id:
                raise AnnotationBuildError(f"Exon at GTF line {line_number} has no gene_id")
            gene_name = attributes.get("gene_name") or None
            record = genes.get(gene_id)
            if record is None:
                genes[gene_id] = GeneRecords(chromosome, gene_name, [(start, end)])
                continue
            if record.chromosome != chromosome:
                raise AnnotationBuildError(
                    f"Gene {gene_id!r} occurs on multiple chromosomes: {record.chromosome!r} and {chromosome!r}"
                )
            if gene_name is not None and record.gene_name is not None and gene_name != record.gene_name:
                raise AnnotationBuildError(
                    f"Gene {gene_id!r} has inconsistent gene_name values: {record.gene_name!r} and {gene_name!r}"
                )
            if record.gene_name is None:
                record.gene_name = gene_name
            record.intervals.append((start, end))

    if not genes:
        raise AnnotationBuildError(f"No exon records found in GTF: {gtf_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_id", "GeneName", "gene_length_bp"])
        for gene_id in sorted(genes):
            record = genes[gene_id]
            writer.writerow([gene_id, record.gene_name or "NA", merge_interval_length(record.intervals)])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build combined_gene_annotation.tsv from GTF exon records")
    parser.add_argument("--gtf", required=True, type=Path, help="Input GTF annotation")
    parser.add_argument("--output", required=True, type=Path, help="Output combined_gene_annotation.tsv")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        build_gene_annotation(args.gtf, args.output)
    except (AnnotationBuildError, OSError) as error:
        raise SystemExit(f"ERROR: {error}") from error
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
