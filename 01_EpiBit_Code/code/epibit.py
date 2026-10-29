#!/usr/bin/env python3
"""
EpiBit: cytosine-resolved information-theoretic analysis of plant epialleles.

Implements the site-wise Bernoulli framework described in the Peloria manuscript:
- pseudocount-corrected methylation probabilities
- binary Shannon entropy
- directional Kullback-Leibler divergence
- Jensen-Shannon divergence
- mean bidirectional KL divergence (Dsym)
- context/locus summaries
- ranked cumulative divergence dimensionality
- pseudocount sensitivity analysis

Python 3.10+
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

REQUIRED_COLUMNS = [
    "site_id", "fragment", "position", "context",
    "WT_methylated", "WT_total", "Peloria_methylated", "Peloria_total"
]

def _log2(x: float) -> float:
    return math.log(x, 2)

def entropy_binary(p: float) -> float:
    """Binary Shannon entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * _log2(p) + (1.0 - p) * _log2(1.0 - p))

def kl_binary(q: float, p: float) -> float:
    """Bernoulli KL divergence D_KL(q || p), in bits."""
    if not (0.0 < p < 1.0 and 0.0 < q < 1.0):
        raise ValueError("KL inputs must be strictly between 0 and 1.")
    return q * _log2(q / p) + (1.0 - q) * _log2((1.0 - q) / (1.0 - p))

def jsd_binary(p: float, q: float) -> float:
    """Jensen-Shannon divergence between Bernoulli(p) and Bernoulli(q), in bits."""
    r = 0.5 * (p + q)
    return 0.5 * kl_binary(p, r) + 0.5 * kl_binary(q, r)

def corrected_probability(methylated: int, total: int, alpha: float = 0.5) -> float:
    """Symmetric pseudocount estimate: (m + alpha) / (n + 2 alpha)."""
    if total < 1:
        raise ValueError("Total informative reads must be >= 1.")
    if methylated < 0 or methylated > total:
        raise ValueError("Methylated reads must satisfy 0 <= m <= n.")
    if alpha <= 0:
        raise ValueError("alpha must be > 0.")
    return (methylated + alpha) / (total + 2.0 * alpha)

def load_matrix(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return df.copy()

def analyze_matrix(df: pd.DataFrame, alpha: float = 0.5) -> pd.DataFrame:
    out = df.copy()
    out["WT_p"] = [
        corrected_probability(int(m), int(n), alpha)
        for m, n in zip(out["WT_methylated"], out["WT_total"])
    ]
    out["Peloria_q"] = [
        corrected_probability(int(m), int(n), alpha)
        for m, n in zip(out["Peloria_methylated"], out["Peloria_total"])
    ]
    out["H_WT_bits"] = out["WT_p"].map(entropy_binary)
    out["H_Peloria_bits"] = out["Peloria_q"].map(entropy_binary)
    out["KL_Peloria_WT_bits"] = [
        kl_binary(q, p) for p, q in zip(out["WT_p"], out["Peloria_q"])
    ]
    out["KL_WT_Peloria_bits"] = [
        kl_binary(p, q) for p, q in zip(out["WT_p"], out["Peloria_q"])
    ]
    out["Dsym_bits"] = 0.5 * (
        out["KL_Peloria_WT_bits"] + out["KL_WT_Peloria_bits"]
    )
    out["JSD_bits"] = [
        jsd_binary(p, q) for p, q in zip(out["WT_p"], out["Peloria_q"])
    ]
    total_dsym = float(out["Dsym_bits"].sum())
    out["Dsym_fraction"] = out["Dsym_bits"] / total_dsym if total_dsym > 0 else 0.0
    return out

def summarize_by_context(site_metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for context, g in site_metrics.groupby("context", dropna=False):
        rows.append({
            "context": context,
            "N": len(g),
            "WT_mean_p": g["WT_p"].mean(),
            "Peloria_mean_q": g["Peloria_q"].mean(),
            "H_WT_bits": g["H_WT_bits"].sum(),
            "H_Peloria_bits": g["H_Peloria_bits"].sum(),
            "KL_Peloria_WT_bits": g["KL_Peloria_WT_bits"].sum(),
            "KL_WT_Peloria_bits": g["KL_WT_Peloria_bits"].sum(),
            "Dsym_bits": g["Dsym_bits"].sum(),
            "JSD_bits": g["JSD_bits"].sum(),
        })
    summary = pd.DataFrame(rows)
    total_kl = summary["KL_Peloria_WT_bits"].sum()
    summary["KL_percent"] = (
        100.0 * summary["KL_Peloria_WT_bits"] / total_kl
        if total_kl > 0 else 0.0
    )
    return summary

def summarize_region(site_metrics: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([{
        "N_cytosines": len(site_metrics),
        "WT_mean_p": site_metrics["WT_p"].mean(),
        "Peloria_mean_q": site_metrics["Peloria_q"].mean(),
        "H_WT_bits": site_metrics["H_WT_bits"].sum(),
        "H_Peloria_bits": site_metrics["H_Peloria_bits"].sum(),
        "KL_Peloria_WT_bits": site_metrics["KL_Peloria_WT_bits"].sum(),
        "KL_WT_Peloria_bits": site_metrics["KL_WT_Peloria_bits"].sum(),
        "Dsym_bits": site_metrics["Dsym_bits"].sum(),
        "JSD_bits": site_metrics["JSD_bits"].sum(),
    }])

def cumulative_dimensionality(site_metrics: pd.DataFrame,
                              thresholds: Iterable[float] = (0.50, 0.80, 0.90, 0.95)):
    ranked = site_metrics.sort_values("Dsym_bits", ascending=False).reset_index(drop=True).copy()
    total = ranked["Dsym_bits"].sum()
    ranked["cumulative_Dsym_fraction"] = ranked["Dsym_bits"].cumsum() / total if total > 0 else 0.0

    dims = []
    for tau in thresholds:
        if total <= 0:
            k = 0
        else:
            k = int(np.searchsorted(ranked["cumulative_Dsym_fraction"].to_numpy(), tau, side="left") + 1)
        dims.append({"threshold": tau, "k_sites": k, "N_sites_total": len(ranked)})
    return ranked, pd.DataFrame(dims)

def run_analysis(matrix_path: str | Path, output_dir: str | Path, alpha: float = 0.5):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    df = load_matrix(matrix_path)
    site = analyze_matrix(df, alpha=alpha)
    context = summarize_by_context(site)
    region = summarize_region(site)
    ranked, dimensionality = cumulative_dimensionality(site)

    site.to_csv(output_dir / "site_metrics.csv", index=False)
    context.to_csv(output_dir / "context_summary.csv", index=False)
    region.to_csv(output_dir / "region_summary.csv", index=False)
    ranked.to_csv(output_dir / "ranked_cytosines.csv", index=False)
    dimensionality.to_csv(output_dir / "dimensionality_thresholds.csv", index=False)

    with pd.ExcelWriter(output_dir / "EpiBit_results.xlsx", engine="openpyxl") as xw:
        site.to_excel(xw, sheet_name="site_metrics", index=False)
        context.to_excel(xw, sheet_name="context_summary", index=False)
        region.to_excel(xw, sheet_name="region_summary", index=False)
        ranked.to_excel(xw, sheet_name="ranked_cytosines", index=False)
        dimensionality.to_excel(xw, sheet_name="dimensionality", index=False)

    return site, context, region, ranked, dimensionality

def sensitivity_analysis(matrix_path: str | Path, output_csv: str | Path,
                         alphas=(0.1, 0.25, 0.5, 1.0, 2.0)):
    df = load_matrix(matrix_path)
    rows = []
    for alpha in alphas:
        site = analyze_matrix(df, alpha=alpha)
        region = summarize_region(site).iloc[0]
        context = summarize_by_context(site).set_index("context")
        _, dims = cumulative_dimensionality(site)
        dim_map = {int(round(100*r.threshold)): int(r.k_sites) for _, r in dims.iterrows()}
        total_kl = float(region["KL_Peloria_WT_bits"])
        chg = float(context.loc["CHG", "KL_Peloria_WT_bits"]) if "CHG" in context.index else 0.0
        chh = float(context.loc["CHH", "KL_Peloria_WT_bits"]) if "CHH" in context.index else 0.0
        rows.append({
            "alpha": alpha,
            "KL_bits": total_kl,
            "Dsym_bits": float(region["Dsym_bits"]),
            "JSD_bits": float(region["JSD_bits"]),
            "CHG_percent_KL": 100.0 * chg / total_kl if total_kl else 0.0,
            "nonCG_percent_KL": 100.0 * (chg + chh) / total_kl if total_kl else 0.0,
            "k50": dim_map.get(50),
            "k80": dim_map.get(80),
            "k90": dim_map.get(90),
            "k95": dim_map.get(95),
        })
    out = pd.DataFrame(rows)
    out.to_csv(output_csv, index=False)
    return out

def main():
    ap = argparse.ArgumentParser(description="EpiBit cytosine-level epiallele analysis")
    sub = ap.add_subparsers(dest="command", required=True)

    a = sub.add_parser("analyze", help="Analyze one cytosine-level methylation matrix")
    a.add_argument("--matrix", required=True)
    a.add_argument("--outdir", required=True)
    a.add_argument("--alpha", type=float, default=0.5)

    s = sub.add_parser("sensitivity", help="Run pseudocount sensitivity analysis")
    s.add_argument("--matrix", required=True)
    s.add_argument("--output", required=True)

    args = ap.parse_args()
    if args.command == "analyze":
        run_analysis(args.matrix, args.outdir, args.alpha)
    elif args.command == "sensitivity":
        sensitivity_analysis(args.matrix, args.output)

if __name__ == "__main__":
    main()
