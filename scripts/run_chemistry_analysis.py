#!/usr/bin/env python3
"""Create descriptive substrate-chemistry tables from analytical replicates."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "chemistry" / "combined_replicate_level_source_data.csv"
OUTPUT = ROOT / "results" / "reproduced" / "chemistry"


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(INPUT)
    summary = (
        source.groupby(["treatment", "variable"], dropna=False).value
        .agg(["count", "mean", "std", "min", "max"])
        .reset_index()
    )
    summary.to_csv(OUTPUT / "descriptive_analytical_summary_long.csv", index=False)

    order = [
        "Nt (%)", "TOC (%)", "TP (ppm)", "C:N", "pH",
        "Total organic matter (%)", "Soluble C (%)", "Hemicellulose (%)",
        "Cellulose (%)", "Lignin (%)",
    ]
    names = {
        "S+": "Full sunlight (S+)",
        "UV-": "Reduced UV (UV-)",
        "S-": "Reduced UV and blue-green light (S-)",
    }
    summary["display"] = summary.apply(
        lambda row: (
            f"{row['mean']:.2f} ± {row['std']:.2f} "
            f"(analytical n={int(row['count'])})"
        ),
        axis=1,
    )
    wide = summary.pivot(index="treatment", columns="variable", values="display")
    wide = wide.reindex(["S+", "UV-", "S-"]).reindex(columns=order)
    wide.insert(0, "Substrate", [names[value] for value in wide.index])
    wide.reset_index(drop=True).to_csv(
        OUTPUT / "Table_1_descriptive_no_inference.csv", index=False
    )

    (OUTPUT / "INTERPRETATION_LIMIT.md").write_text(
        "# Interpretation limit\n\n"
        "Each sunlight treatment was represented by one pooled exposed-litter batch. "
        "Repeated measurements are analytical replicates or aliquots, not independent "
        "exposure-level biological replicates. Table 1 therefore reports descriptive "
        "means and analytical variation only; treatment ANOVA and causal treatment "
        "claims are not supported by this design.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
