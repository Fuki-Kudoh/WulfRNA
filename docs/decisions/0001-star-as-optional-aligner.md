# ADR 0001: Add STAR as an optional alignment branch

## Status

Accepted

## Context

WulfRNA currently treats Salmon and kallisto as transcript
quantification backends. STAR produces genome alignments and gene-level
counts but does not provide an equivalent transcript TPM output.

## Options considered

1. Add STAR to --quantifier.
2. Replace Salmon with STAR.
3. Add STAR as an independent optional aligner.

## Decision

Use an independent optional alignment phase.

## Rationale

This preserves the existing quantification contract while adding BAM
and splice-junction outputs for inspection and downstream analysis.

## Consequences

- The phase list becomes conditional.
- Alignment and quantification have separate resume semantics.
- STAR counts are not combined with Salmon/kallisto matrices initially.
