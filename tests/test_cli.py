import pytest

from wulfrna.cli import parse_args


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


def test_aligner_defaults_to_none_preserving_legacy_invocation():
    args = parse_args([
        "work",
        "--reference",
        "ref",
        "--stranded",
        "reverse",
        "--threads",
        "4",
    ])

    assert args.command == "run"
    assert args.quantifier == "salmon"
    assert args.aligner == "none"


def test_aligner_accepts_star():
    args = parse_args([*base_args(), "--aligner", "star"])

    assert args.aligner == "star"
    assert args.quantifier == "salmon"


def test_aligner_rejects_unknown_value():
    with pytest.raises(SystemExit) as exc:
        parse_args([*base_args(), "--aligner", "bowtie2"])

    assert exc.value.code == 2
