#!/usr/bin/env python3

"""Apply CAP's conservative source-preservation policy to review candidates."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

RATIONALES = {
    "gpl-source-review": (
        "yes",
        "Distributed GPL-family runtime component; preserve exact corresponding "
        "source, recipe, patches, build scripts, and notices.",
    ),
    "lgpl-source-review": (
        "yes",
        "Conservative LGPL compliance: preserve exact source and build material "
        "regardless of dynamic-link and relinking analysis.",
    ),
    "alternative-license-review": (
        "yes",
        "A permissive or weaker-copyleft alternative may apply, but preserve exact "
        "source and notices conservatively rather than relying on omission.",
    ),
    "exception-review": (
        "yes",
        "A runtime/linking exception may permit distribution without source, but "
        "preserve exact source and notices conservatively.",
    ),
}

EXPECTED = {
    "alternative-license-review": {
        "cairo", "gmp", "libev", "libfreetype", "libfreetype6", "perl",
    },
    "exception-review": {
        "libgcc", "libgcc-ng", "libgfortran", "libgfortran-ng",
        "libgfortran5", "libgomp", "libstdcxx", "libstdcxx-ng", "openjdk",
    },
    "gpl-source-review": {
        "bioconductor-msa", "coreutils", "gawk", "glpk", "gsl",
        "ld_impl_linux-64", "libgettextpo", "r-ade4", "r-ape", "r-base",
        "r-codetools", "r-digest", "r-doparallel", "r-getopt", "r-lattice",
        "r-mass", "r-mime", "r-nlme", "r-pixmap", "r-rcpp", "r-rlang",
        "r-segmented", "r-seqinr", "r-sp", "r-stringdist", "readline", "sed",
    },
    "lgpl-source-review": {
        "alsa-lib", "fribidi", "graphite2", "keyutils", "libasprintf",
        "libglib", "libiconv", "libnsl", "libxcrypt", "mpfr", "pango",
    },
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    with args.input.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = reader.fieldnames
    if not rows or fields is None:
        raise SystemExit("Review CSV is empty")

    actual: dict[str, set[str]] = {}
    for row in rows:
        review_class = row["review_class"]
        actual.setdefault(review_class, set()).add(row["name"])
    if actual != EXPECTED:
        missing = {
            key: sorted(names - actual.get(key, set()))
            for key, names in EXPECTED.items() if names - actual.get(key, set())
        }
        unexpected = {
            key: sorted(names - EXPECTED.get(key, set()))
            for key, names in actual.items() if names - EXPECTED.get(key, set())
        }
        raise SystemExit(
            f"Candidate set changed; review policy before proceeding. "
            f"Missing={missing}, unexpected={unexpected}"
        )

    for row in rows:
        decision, rationale = RATIONALES[row["review_class"]]
        row["collect_source"] = decision
        row["review_note"] = rationale

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
