#!/usr/bin/env python3
"""Publish the self-contained customer-facing report."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        shutil.rmtree(args.output)
    args.output.mkdir(parents=True, exist_ok=True)
    source = args.scratch / "report/final_ISO_SLAM_AS_report.html"
    if not source.is_file():
        raise FileNotFoundError(source)
    shutil.copy2(source, args.output / source.name)
    (args.output / "stage_complete.txt").write_text("AS report complete\n", encoding="utf-8")


if __name__ == "__main__":
    main()
