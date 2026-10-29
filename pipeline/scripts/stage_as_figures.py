#!/usr/bin/env python3
"""Publish figures produced by the AS core calculation."""
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
    source = args.scratch / "figures"
    if not source.is_dir():
        raise FileNotFoundError(source)
    shutil.copytree(source, args.output / "figures")
    (args.output / "stage_complete.txt").write_text("AS figures complete\n", encoding="utf-8")


if __name__ == "__main__":
    main()
