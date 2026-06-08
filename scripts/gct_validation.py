from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr, spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
EXTERNAL_DIR = PROJECT_ROOT / "data" / "external" / "gct"
FIGURES_DIR = PROJECT_ROOT / "figures"
GCT_SPECIES = EXTERNAL_DIR / "gct_species.csv"


def load_gct_species(path=GCT_SPECIES):
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Run `python scripts/import_gct_dataset.py` first."
        )

    df = pd.read_csv(path)
    required = {"species", "GCw", "Topt", "thermal_group"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {path}: {sorted(missing)}")

    return df.dropna(subset=["GCw", "Topt"]).reset_index(drop=True)


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    df = load_gct_species()

    pearson_r, pearson_p = pearsonr(df["GCw"], df["Topt"])
    spearman_r, spearman_p = spearmanr(df["GCw"], df["Topt"])

    print("=== gcT External Validation: GCw vs Topt ===")
    print(f"Rows: {len(df)}")
    print(f"Pearson r : {pearson_r:.4f} (p={pearson_p:.2e})")
    print(f"Spearman r: {spearman_r:.4f} (p={spearman_p:.2e})")

    direction = "positive" if pearson_r > 0 else "negative"
    print(
        "Interpretation: gcT shows a "
        f"{direction} genomic GC/Topt association in this external dataset."
    )

    plt.figure(figsize=(9, 6))
    sns.scatterplot(
        data=df,
        x="GCw",
        y="Topt",
        hue="thermal_group",
        alpha=0.65,
        s=35,
        edgecolor="white",
        linewidth=0.3,
    )
    sns.regplot(
        data=df,
        x="GCw",
        y="Topt",
        scatter=False,
        color="black",
        line_kws={"linestyle": "--", "linewidth": 1.5},
    )
    plt.title("gcT External Validation: Genomic GC Content vs Topt", fontweight="bold")
    plt.xlabel("Genomic GC content (%)")
    plt.ylabel("Topt (degC)")
    plt.tight_layout()
    out_path = FIGURES_DIR / "gct_external_gc_topt.png"
    plt.savefig(out_path, bbox_inches="tight")
    plt.close()
    print(f"Saved figure: {out_path}")


if __name__ == "__main__":
    main()
