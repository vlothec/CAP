#!/usr/bin/env python3

"""Apply reviewed resolutions to a Conda source-notice manifest."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

KEYS = ("name", "version", "build")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source_notice_manifest", type=Path)
    parser.add_argument("review_policy", type=Path)
    parser.add_argument("shared_license_directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Output path already exists: {args.output}")

    rows = read_csv(args.source_notice_manifest)
    policies = read_csv(args.review_policy)
    policy_by_key = {tuple(row[key] for key in KEYS): row for row in policies}
    if len(policy_by_key) != len(policies):
        raise SystemExit("Duplicate package key in review policy")

    resolved: list[dict[str, str]] = []
    used: set[tuple[str, str, str]] = set()
    for row in rows:
        key = tuple(row[field] for field in KEYS)
        status = row["review_status"]
        method = status
        evidence = row["source_license_files"]
        note = "License-like files extracted from verified source archive."
        if status == "REVIEW_REQUIRED":
            policy = policy_by_key.get(key)
            if policy is None:
                raise SystemExit(f"No reviewed resolution for {'-'.join(key)}")
            used.add(key)
            method = policy["resolution"]
            evidence = policy["canonical_license_files"]
            note = policy["review_note"]
            if method == "R-shared-license":
                for filename in evidence.split(" | "):
                    if not (args.shared_license_directory / filename).is_file():
                        raise SystemExit(f"Missing canonical R license: {filename}")
            elif method not in {
                "recipe-license-declaration",
            }:
                raise SystemExit(f"Unknown reviewed resolution: {method}")
        resolved.append({
            **row,
            "final_resolution": method,
            "final_notice_evidence": evidence,
            "final_review_note": note,
        })

    unused = set(policy_by_key) - used
    if unused:
        raise SystemExit(
            "Review policy contains unused package keys: "
            + ", ".join("-".join(key) for key in sorted(unused))
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(resolved[0])
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(resolved)
    print(f"Resolved source-level notice review for {len(resolved)} packages")


if __name__ == "__main__":
    main()
