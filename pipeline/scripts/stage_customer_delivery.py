#!/usr/bin/env python3
"""Create the strict customer delivery from published AS outputs."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path


FIGURES = [
    "figures/slam_qc/new_rna_fraction.png", "figures/slam_qc/tc_conversion_qc.png",
    "figures/total_rna/highres/event_candidate_counts.png", "figures/total_rna/highres/delta_psi_heatmap.png",
    "figures/new_rna/new_rna_event_coverage.png", "figures/new_rna/new_rna_as_candidate_counts.png",
    "figures/new_rna/new_rna_event_type_candidate_counts.png", "figures/new_rna/new_rna_delta_psi_heatmap.png",
    "figures/integration/total_vs_new_delta_psi_6panel.png", "figures/integration/total_new_delta_psi_heatmap.png",
    "figures/integration/event_type_comparison.png",
] + [f"figures/case_studies/case_{i:02d}.png" for i in range(1, 7)]
TABLES = [
    "tables/case_studies.tsv", "tables/new_rna/as_event_summary.tsv", "tables/new_rna/event_psi.tsv",
    "tables/new_rna/event_support.tsv", "tables/total_rna/event_psi.tsv",
] + [f"tables/new_rna/contrasts/{name}.event_candidates.tsv" for name in ("A_vs_NC", "B_vs_NC", "AB_vs_NC", "A_vs_B", "AB_vs_A", "AB_vs_B")] + [f"tables/integration/{name}.total_vs_new.tsv" for name in ("A_vs_NC", "B_vs_NC", "AB_vs_NC", "A_vs_B", "AB_vs_A", "AB_vs_B")]
METHODS = ["contrast_direction.md", "event_definition.tsv", "new_rna_psi_filtering.md", "sample_metadata.tsv"]


def copy_one(source: Path, target: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        shutil.rmtree(args.output)
    args.output.mkdir(parents=True, exist_ok=True)
    for rel in FIGURES:
        copy_one(args.figures / rel, args.output / rel)
    for rel in TABLES:
        copy_one(args.tables / rel, args.output / rel)
    for rel in METHODS:
        copy_one(args.tables / "provenance" / rel, args.output / "methods" / rel)
    copy_one(args.report / "final_ISO_SLAM_AS_report.html", args.output / "report/final_ISO_SLAM_AS_report.html")
    (args.output / "README.md").write_text("# ISO-SLAM-seq 可变剪切分析结果\n\n本目录提供中文离线分析报告、报告图表、核心可变剪切结果表和方法说明。样本标签为匿名分析编号。\n\n- 报告：`report/final_ISO_SLAM_AS_report.html`\n- 图表：`figures/`\n- 核心结果：`tables/`\n- 方法和样本说明：`methods/`\n", encoding="utf-8")


if __name__ == "__main__":
    main()
