"""Design-adjusted, all-genus differential-abundance analysis.

Top-12 abundance is used only for composition figures. Every genus with at least 10
counts across the retained samples is tested here. PyDESeq2 implements the DESeq2 model;
blocking variables are nuisance fixed effects because negative-binomial DESeq2 models do
not support random effects.
"""
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

REPO_ROOT = Path(__file__).resolve().parents[1]
ROOT = REPO_ROOT / "results" / "reproduced" / "amplicons"
COLORS = {"S+":"#E3AC18","UV-":"#77BCE8","S-":"#555555",
          "SusS+":"#E3AC18","SusUV-":"#77BCE8","SusS-":"#555555"}

def analyze(experiment: str, marker: str):
    counts = pd.read_csv(ROOT / f"{experiment}_{marker}_all_genus_counts.csv").set_index("sample_id")
    meta = pd.read_csv(ROOT / f"{experiment}_{marker}_metadata.csv").set_index("sample_id")
    common = counts.index.intersection(meta.index)
    counts = counts.loc[common].round().astype(int)
    meta = meta.loc[common].copy()
    keep = counts.sum(axis=0) >= 10
    counts = counts.loc[:, keep]
    meta["Treatment"] = meta["Treatment"].astype(str)
    if experiment == "field":
        levels = ["S+","UV-","S-"]
        meta["Block"] = meta["Block"].astype(str)
        design = "~Block + Treatment"
    else:
        levels = ["SusS+","SusUV-","SusS-"]
        # Inoculum-by-replicate block subsumes both nuisance variables.
        meta["PermutationBlock"] = meta["PermutationBlock"].astype(str)
        design = "~PermutationBlock + Treatment"
    meta["Treatment"] = pd.Categorical(meta["Treatment"], categories=levels)
    dds = DeseqDataSet(counts=counts, metadata=meta, design=design,
                       refit_cooks=True, n_cpus=1, quiet=True)
    dds.deseq2()
    comparisons = [(levels[0],levels[1]),(levels[0],levels[2]),(levels[1],levels[2])]
    all_results=[]
    for numerator, denominator in comparisons:
        stat = DeseqStats(dds, contrast=["Treatment", numerator, denominator],
                          alpha=.05, independent_filter=True, quiet=True, n_cpus=1)
        stat.summary()
        res = stat.results_df.copy().reset_index().rename(columns={"index":"genus"})
        res["experiment"] = experiment; res["marker"] = marker
        res["contrast"] = f"{numerator}_vs_{denominator}"
        res["numerator"] = numerator; res["denominator"] = denominator
        res["total_counts"] = res["genus"].map(counts.sum(axis=0))
        all_results.append(res)
    out=pd.concat(all_results,ignore_index=True)
    out.to_csv(ROOT/f"{experiment}_{marker}_all_genus_deseq2.csv",index=False)
    return out, len(counts), counts.shape[1], design

summary=[]; results={}
for exp in ["field","microcosm"]:
    for marker in ["ITS","16S"]:
        out,n,ng,design=analyze(exp,marker)
        results[(exp,marker)]=out
        summary.append({"experiment":exp,"marker":marker,"N":n,"genera_tested":ng,
                        "design":design,"significant_rows_q_lt_0.05":int((out.padj<.05).sum())})
pd.DataFrame(summary).to_csv(ROOT/"amplicon_deseq2_summary.csv",index=False)
combined=pd.concat(results.values(),ignore_index=True)
combined.loc[combined.padj.lt(.05).fillna(False)].sort_values(
    ["experiment","marker","contrast","padj","genus"]
).to_csv(ROOT/"all_experiments_significant_genera.csv",index=False)

def readable_taxon(value: str) -> str:
    value = value.replace("_gen_Incertae_sedis", " genus incertae sedis")
    value = value.replace("_", " ")
    return "\n".join(textwrap.wrap(value, width=34, break_long_words=False))


def draw_volcano(ax, d: pd.DataFrame, marker: str, panel: str, contrast: str,
                 label_limit: int | None = None) -> None:
    d = d.copy()
    d["minuslog10q"] = -np.log10(d.padj.clip(lower=np.nextafter(0, 1)))
    sig = d.padj.lt(.05).fillna(False)
    ax.scatter(d.loc[~sig, "log2FoldChange"], d.loc[~sig, "minuslog10q"],
               s=30, c="#AFAFAF", edgecolor="black", linewidth=.35, alpha=.68)
    ax.scatter(d.loc[sig, "log2FoldChange"], d.loc[sig, "minuslog10q"],
               s=52, c="#B3261E", edgecolor="black", linewidth=.5, alpha=.92, zorder=3)
    ax.axhline(-np.log10(.05), ls="--", lw=.85, c="#333333")
    ax.axvline(0, lw=.65, c="#333333")

    # Only significant categories are labelled. Text is placed immediately
    # beside its point; no connector lines are used.
    lab = d.loc[sig].sort_values(["padj", "total_counts"], ascending=[True, False])
    if label_limit is not None:
        lab = lab.head(label_limit)
    offsets = [(7, 7), (7, -13), (7, 8), (7, -12), (-7, 8), (-7, -12)]
    field_offsets = {
        # Separate the three bacterial labels clustered near q = 0.003 without
        # moving them away from their points or adding connector lines.
        "Modestobacter": (8, 11, "left"),
        "Psychroglaciecola": (-8, -7, "right"),
        "Naasia": (8, -18, "left"),
        "1174-901-12": (8, -12, "left"),
        "Microvirga": (-8, 8, "right"),
        "Dothideomycetes_gen_Incertae_sedis": (8, 8, "left"),
    }
    for i, (_, row) in enumerate(lab.iterrows()):
        if str(row.genus) in field_offsets:
            dx, dy, ha = field_offsets[str(row.genus)]
        else:
            dx, dy = offsets[i % len(offsets)]
            ha = "left" if row.log2FoldChange >= 0 else "right"
            dx = abs(dx) if ha == "left" else -abs(dx)
        ax.annotate(readable_taxon(str(row.genus)),
                    (row.log2FoldChange, row.minuslog10q),
                    xytext=(dx, dy), textcoords="offset points", fontsize=10.5,
                    ha=ha, va="bottom" if dy >= 0 else "top",
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": .72, "pad": .6})

    ax.set_title(f"{panel}  {'Fungi' if marker == 'ITS' else 'Bacteria'}", loc="left",
                 fontsize=15, fontweight="bold")
    ax.set_xlabel("log$_2$ fold change (" + contrast.replace("_vs_", " / ") + ")", fontsize=12)
    ax.set_ylabel("−log$_{10}$(BH-adjusted p)", fontsize=12)
    ax.tick_params(axis="both", labelsize=10.5)
    ax.spines[["top", "right"]].set_visible(False)


for exp in ["field", "microcosm"]:
    if exp == "field":
        # Publication Figure S2: all three field contrasts and every eligible
        # category. Only significant results are labelled.
        levels = ["S+_vs_UV-", "S+_vs_S-", "UV-_vs_S-"]
        fig, axes = plt.subplots(2, 3, figsize=(19, 10.5), sharex=False, sharey=False)
        for r, marker in enumerate(["ITS", "16S"]):
            df = results[(exp, marker)]
            for c, contrast in enumerate(levels):
                panel = "ABC"[c] if r == 0 else "DEF"[c]
                ax = axes[r, c]
                draw_volcano(ax, df.loc[df.contrast.eq(contrast)], marker,
                             panel, contrast)
                comparison = contrast.replace("_vs_", " vs ").replace("UV-", "UV−").replace("S-", "S−")
                community = "Fungi" if marker == "ITS" else "Bacteria"
                ax.set_title(f"{panel}  {community}: {comparison}", loc="left",
                             fontsize=14, fontweight="bold")
        fig.suptitle("Field sunlight exposure: genus-level differential abundance",
                     fontsize=18, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, .96))
    else:
        levels = ["SusS+_vs_SusUV-", "SusS+_vs_SusS-", "SusUV-_vs_SusS-"]
        fig, axes = plt.subplots(2, 3, figsize=(18, 10.5), sharex=False, sharey=False)
        for r, marker in enumerate(["ITS", "16S"]):
            df = results[(exp, marker)]
            for c, contrast in enumerate(levels):
                draw_volcano(axes[r, c], df.loc[df.contrast.eq(contrast)], marker,
                             "ABC"[c] if r == 0 else "DEF"[c], contrast,
                             label_limit=6)
        fig.suptitle("Microcosm litter photodegradation: all-genus differential abundance",
                     fontsize=18, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, .96))
    fig.savefig(ROOT / f"{exp}_all_genus_volcano.pdf", bbox_inches="tight")
    fig.savefig(ROOT / f"{exp}_all_genus_volcano.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
