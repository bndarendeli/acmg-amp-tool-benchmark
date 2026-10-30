"""Build the variant-level evidence concordance summaries used in Figure 6."""

from pathlib import Path

import pandas as pd


OUTPUT_DIR = Path(__file__).resolve().parent / "outputs"
MASTER_FILE = OUTPUT_DIR / "fig6_master_variant_tool.csv"
OUTPUT_FILE = OUTPUT_DIR / "fig6_evidence_distributions.csv"


def summarize(dataset: str, tool: str, rows: pd.DataFrame) -> dict:
    jaccard = rows["criterion_jaccard"].dropna().astype(float)
    n = len(jaccard)
    if n == 0:
        raise ValueError(f"No evidence-evaluable rows for {dataset}/{tool}")
    return {
        "dataset": dataset,
        "tool": tool,
        "n_classification_concordant_with_evaluable_evidence": n,
        "exact_evidence_match_rate": rows.loc[jaccard.index, "exact_evidence_match"].astype(bool).mean(),
        "mean_criterion_jaccard": jaccard.mean(),
        "median_criterion_jaccard": jaccard.median(),
        "q1_criterion_jaccard": jaccard.quantile(0.25),
        "q3_criterion_jaccard": jaccard.quantile(0.75),
        "jaccard_1.00": (jaccard == 1.0).mean(),
        "jaccard_0.75_to_1.00": ((jaccard >= 0.75) & (jaccard < 1.0)).mean(),
        "jaccard_0.50_to_0.75": ((jaccard >= 0.50) & (jaccard < 0.75)).mean(),
        "jaccard_0.25_to_0.50": ((jaccard >= 0.25) & (jaccard < 0.50)).mean(),
        "jaccard_0.00_to_0.25": ((jaccard > 0.00) & (jaccard < 0.25)).mean(),
        "jaccard_0.00": (jaccard == 0.0).mean(),
    }


def main() -> None:
    master = pd.read_csv(MASTER_FILE, low_memory=False)
    eligible = master.loc[
        master["classification_match"].eq(True) & master["criterion_jaccard"].notna()
    ].copy()
    records = []
    for dataset, dataset_rows in eligible.groupby("dataset", sort=False):
        for tool, tool_rows in dataset_rows.groupby("tool", sort=True):
            records.append(summarize(dataset, tool, tool_rows))
        records.append(summarize(dataset, "POOLED", dataset_rows))

    result = pd.DataFrame(records)
    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(OUTPUT_FILE, index=False)
    print(f"Saved {len(result)} dataset/tool summaries: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
