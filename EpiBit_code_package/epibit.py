#!/usr/bin/env python3
"""
EpiBit: cytosine-level information-theoretic analysis of plant epialleles
from targeted bisulfite/NGS reads or from a precomputed cytosine count matrix.

Reference implementation used to document the manuscript analysis.
Language: Python 3.10+
Dependencies: numpy, pandas, openpyxl (only for Excel export)

Main modes
----------
1) from-matrix: compute EpiBit metrics from a cytosine-level count matrix.
2) align:       parse bisulfite FASTA reads, align to reference fragments, call
                methylation counts, and compute EpiBit metrics.

Important analytical defaults for the manuscript implementation
---------------------------------------------------------------
MIN_IDENTITY = 0.85
ALPHA        = 0.5
MIN_COVERAGE = None  # no additional post-alignment coverage threshold

The alignment score allows the expected bisulfite C->T pattern: at a cytosine
in the reference, both C and T in the read are treated as compatible with the
reference for alignment scoring. For methylation calling, C is counted as
methylated and T as unmethylated. Other bases at reference-C positions are not
included in methylated/total counts.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

MIN_IDENTITY = 0.85
ALPHA = 0.5
MIN_COVERAGE = None
DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path: str | Path) -> List[Tuple[str, str]]:
    """Read a FASTA file and return [(record_id, sequence), ...]."""
    records: List[Tuple[str, str]] = []
    current_id: Optional[str] = None
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records.append((current_id, "".join(chunks).upper()))
                current_id = line[1:].strip() or f"record_{len(records)+1}"
                chunks = []
            else:
                chunks.append(line)
    if current_id is not None:
        records.append((current_id, "".join(chunks).upper()))
    return records


def reverse_complement(seq: str) -> str:
    return seq.translate(DNA_COMPLEMENT)[::-1].upper()


def compatible_for_bisulfite_alignment(ref_base: str, read_base: str) -> bool:
    """Return True if read_base is compatible with ref_base under bisulfite logic."""
    ref_base = ref_base.upper()
    read_base = read_base.upper()
    if read_base == "N":
        return False
    if ref_base == "C":
        return read_base in {"C", "T"}
    return read_base == ref_base


@dataclass
class AlignmentResult:
    accepted: bool
    orientation: str
    offset: int
    identity: float
    aligned_length: int
    matches: int
    mismatches: int
    oriented_sequence: str


def score_alignment(read_seq: str, ref_seq: str, offset: int) -> Tuple[float, int, int, int]:
    """
    Score read_seq against ref_seq at a given offset.

    offset is the reference coordinate corresponding to read position 0.
    Negative offsets allow the read to start before the reference.
    """
    matches = 0
    mismatches = 0
    compared = 0
    ref_len = len(ref_seq)
    for j, rb in enumerate(read_seq):
        k = offset + j
        if k < 0 or k >= ref_len:
            continue
        if rb.upper() == "N":
            continue
        compared += 1
        if compatible_for_bisulfite_alignment(ref_seq[k], rb):
            matches += 1
        else:
            mismatches += 1
    identity = matches / compared if compared else 0.0
    return identity, compared, matches, mismatches


def best_bisulfite_alignment(
    read_seq: str,
    ref_seq: str,
    min_identity: float = MIN_IDENTITY,
    min_overlap: int = 50,
    allow_reverse_complement: bool = True,
) -> AlignmentResult:
    """Find the best direct/reverse-complement local placement of a read."""
    candidates = [("direct", read_seq.upper())]
    if allow_reverse_complement:
        candidates.append(("rc", reverse_complement(read_seq)))

    best: Optional[AlignmentResult] = None
    # scan all placements with at least a possible overlap of min_overlap
    for orientation, seq in candidates:
        start = -len(seq) + min_overlap
        end = len(ref_seq) - min_overlap
        for offset in range(start, end + 1):
            identity, aligned_len, matches, mismatches = score_alignment(seq, ref_seq, offset)
            if aligned_len < min_overlap:
                continue
            result = AlignmentResult(
                accepted=identity >= min_identity,
                orientation=orientation,
                offset=offset,
                identity=identity,
                aligned_length=aligned_len,
                matches=matches,
                mismatches=mismatches,
                oriented_sequence=seq,
            )
            if best is None:
                best = result
            else:
                key = (result.identity, result.aligned_length, result.matches)
                best_key = (best.identity, best.aligned_length, best.matches)
                if key > best_key:
                    best = result
    if best is None:
        return AlignmentResult(False, "none", 0, 0.0, 0, 0, 0, read_seq.upper())
    best.accepted = best.identity >= min_identity
    return best


def cytosine_context(ref_seq: str, pos0: int) -> str:
    """Classify reference cytosine at 0-based pos0 as CG, CHG or CHH."""
    seq = ref_seq.upper()
    if seq[pos0] != "C":
        raise ValueError("cytosine_context called on a non-C position")
    if pos0 + 1 < len(seq) and seq[pos0 + 1] == "G":
        return "CG"
    if pos0 + 2 < len(seq) and seq[pos0 + 2] == "G" and seq[pos0 + 1] != "G":
        return "CHG"
    return "CHH"


def annotate_cytosines(ref_seq: str, fragment: str, locus: str = "") -> pd.DataFrame:
    rows = []
    for pos0, base in enumerate(ref_seq.upper()):
        if base == "C":
            rows.append({
                "site_id": f"{locus + '_' if locus else ''}{fragment}_C{len(rows)+1:03d}",
                "fragment": fragment,
                "position": pos0 + 1,
                "context": cytosine_context(ref_seq, pos0),
            })
    return pd.DataFrame(rows)


def call_counts_for_sample(
    records: List[Tuple[str, str]],
    ref_seq: str,
    annotation: pd.DataFrame,
    min_identity: float = MIN_IDENTITY,
    min_overlap: int = 50,
    allow_reverse_complement: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Align reads and count methylated/total observations at each annotated cytosine."""
    counts = annotation[["site_id", "fragment", "position", "context"]].copy()
    counts["methylated_reads"] = 0
    counts["total_reads"] = 0
    counts["ignored_non_CT_at_C_sites"] = 0
    pos_to_idx = {int(row.position) - 1: idx for idx, row in counts.iterrows()}
    qc_rows = []

    for read_id, seq in records:
        aln = best_bisulfite_alignment(
            seq, ref_seq,
            min_identity=min_identity,
            min_overlap=min_overlap,
            allow_reverse_complement=allow_reverse_complement,
        )
        qc_rows.append({
            "read_id": read_id,
            "accepted": aln.accepted,
            "orientation": aln.orientation,
            "offset": aln.offset,
            "identity": aln.identity,
            "aligned_length": aln.aligned_length,
            "matches": aln.matches,
            "mismatches": aln.mismatches,
        })
        if not aln.accepted:
            continue
        seq_oriented = aln.oriented_sequence
        for pos0, idx in pos_to_idx.items():
            j = pos0 - aln.offset
            if j < 0 or j >= len(seq_oriented):
                continue
            rb = seq_oriented[j].upper()
            if rb == "C":
                counts.at[idx, "methylated_reads"] += 1
                counts.at[idx, "total_reads"] += 1
            elif rb == "T":
                counts.at[idx, "total_reads"] += 1
            elif rb != "N":
                counts.at[idx, "ignored_non_CT_at_C_sites"] += 1
    return counts, pd.DataFrame(qc_rows)


def binary_entropy(x: float) -> float:
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log2(x) - (1.0 - x) * math.log2(1.0 - x)


def bernoulli_kl(a: float, b: float) -> float:
    return a * math.log2(a / b) + (1.0 - a) * math.log2((1.0 - a) / (1.0 - b))


def add_epibit_metrics(matrix: pd.DataFrame, alpha: float = ALPHA) -> pd.DataFrame:
    """Add corrected probabilities and site-wise information metrics to a count matrix."""
    df = matrix.copy()
    wt_m = df["WT_methylated_reads"].astype(float)
    wt_n = df["WT_total_reads"].astype(float)
    pel_m = df["Peloria_methylated_reads"].astype(float)
    pel_n = df["Peloria_total_reads"].astype(float)

    df["WT_p_corrected"] = (wt_m + alpha) / (wt_n + 2.0 * alpha)
    df["Peloria_q_corrected"] = (pel_m + alpha) / (pel_n + 2.0 * alpha)
    df["Delta_q_minus_p"] = df["Peloria_q_corrected"] - df["WT_p_corrected"]
    df["H_WT_bits"] = df["WT_p_corrected"].map(binary_entropy)
    df["H_Peloria_bits"] = df["Peloria_q_corrected"].map(binary_entropy)
    df["D_KL_Peloria_vs_WT_bits"] = [
        bernoulli_kl(q, p) for p, q in zip(df["WT_p_corrected"], df["Peloria_q_corrected"])
    ]
    df["D_KL_WT_vs_Peloria_bits"] = [
        bernoulli_kl(p, q) for p, q in zip(df["WT_p_corrected"], df["Peloria_q_corrected"])
    ]
    df["D_sym_bits"] = 0.5 * (df["D_KL_Peloria_vs_WT_bits"] + df["D_KL_WT_vs_Peloria_bits"])
    midpoint = 0.5 * (df["WT_p_corrected"] + df["Peloria_q_corrected"])
    df["JSD_bits"] = [
        0.5 * bernoulli_kl(q, r) + 0.5 * bernoulli_kl(p, r)
        for p, q, r in zip(df["WT_p_corrected"], df["Peloria_q_corrected"], midpoint)
    ]
    return df


def summarize_by_context(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for context, sub in metrics.groupby("context", sort=False):
        rows.append({
            "context": context,
            "N": len(sub),
            "WT_mean_p": sub["WT_p_corrected"].mean(),
            "Peloria_mean_q": sub["Peloria_q_corrected"].mean(),
            "H_WT_bits": sub["H_WT_bits"].sum(),
            "H_Peloria_bits": sub["H_Peloria_bits"].sum(),
            "KL_Peloria_vs_WT_bits": sub["D_KL_Peloria_vs_WT_bits"].sum(),
            "KL_WT_vs_Peloria_bits": sub["D_KL_WT_vs_Peloria_bits"].sum(),
            "D_sym_bits": sub["D_sym_bits"].sum(),
            "JSD_bits": sub["JSD_bits"].sum(),
        })
    total = {
        "context": "Total",
        "N": len(metrics),
        "WT_mean_p": metrics["WT_p_corrected"].mean(),
        "Peloria_mean_q": metrics["Peloria_q_corrected"].mean(),
        "H_WT_bits": metrics["H_WT_bits"].sum(),
        "H_Peloria_bits": metrics["H_Peloria_bits"].sum(),
        "KL_Peloria_vs_WT_bits": metrics["D_KL_Peloria_vs_WT_bits"].sum(),
        "KL_WT_vs_Peloria_bits": metrics["D_KL_WT_vs_Peloria_bits"].sum(),
        "D_sym_bits": metrics["D_sym_bits"].sum(),
        "JSD_bits": metrics["JSD_bits"].sum(),
    }
    rows.append(total)
    out = pd.DataFrame(rows)
    total_dsym = float(metrics["D_sym_bits"].sum())
    out["D_sym_percent"] = np.where(out["context"] == "Total", 100.0, out["D_sym_bits"] / total_dsym * 100.0)
    return out


def dimensionality(metrics: pd.DataFrame, thresholds=(0.5, 0.8, 0.9, 0.95)) -> Tuple[pd.DataFrame, pd.DataFrame]:
    ranked = metrics.sort_values("D_sym_bits", ascending=False).reset_index(drop=True).copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    total = ranked["D_sym_bits"].sum()
    ranked["relative_Dsym"] = ranked["D_sym_bits"] / total if total > 0 else np.nan
    ranked["cumulative_Dsym_fraction"] = ranked["relative_Dsym"].cumsum()
    rows = []
    for t in thresholds:
        k = int((ranked["cumulative_Dsym_fraction"] >= t).idxmax() + 1) if total > 0 else np.nan
        rows.append({"threshold": f"{int(t*100)}%", "k_sites": k, "N_sites_total": len(ranked), "fraction_sites": k / len(ranked)})
    return ranked, pd.DataFrame(rows)


def apply_min_coverage(df: pd.DataFrame, min_coverage: Optional[int]) -> pd.DataFrame:
    if min_coverage is None:
        return df
    return df[(df["WT_total_reads"] >= min_coverage) & (df["Peloria_total_reads"] >= min_coverage)].copy()


def run_from_matrix(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    matrix = pd.read_csv(args.matrix)
    filtered = apply_min_coverage(matrix, args.min_coverage)
    metrics = add_epibit_metrics(filtered, alpha=args.alpha)
    summary = summarize_by_context(metrics)
    ranked, thresholds = dimensionality(metrics)
    metrics.to_csv(outdir / "derived_metrics.csv", index=False)
    summary.to_csv(outdir / "summary_by_context.csv", index=False)
    ranked.to_csv(outdir / "ranked_dimensionality.csv", index=False)
    thresholds.to_csv(outdir / "dimensionality_thresholds.csv", index=False)
    try:
        with pd.ExcelWriter(outdir / "EpiBit_results.xlsx", engine="openpyxl") as writer:
            filtered.to_excel(writer, sheet_name="cytosine_matrix", index=False)
            metrics.to_excel(writer, sheet_name="derived_metrics", index=False)
            summary.to_excel(writer, sheet_name="summary_by_context", index=False)
            ranked.to_excel(writer, sheet_name="ranked_dimensionality", index=False)
            thresholds.to_excel(writer, sheet_name="thresholds", index=False)
    except Exception:
        pass


def run_align(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    with open(args.config, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    all_matrices = []
    qc_summary_rows = []
    for locus in config["loci"]:
        locus_name = locus["locus"]
        fragment = locus["fragment"]
        ref_seq = locus["reference_sequence"].upper()
        ann = annotate_cytosines(ref_seq, fragment=fragment, locus=locus_name)
        sample_counts: Dict[str, pd.DataFrame] = {}
        for sample in locus["samples"]:
            condition = sample["condition"]
            records = read_fasta(sample["fasta"])
            counts, read_qc = call_counts_for_sample(
                records, ref_seq, ann,
                min_identity=args.min_identity,
                min_overlap=args.min_overlap,
                allow_reverse_complement=sample.get("allow_reverse_complement", True),
            )
            sample_counts[condition] = counts
            read_qc.to_csv(outdir / f"{locus_name}_{fragment}_{condition}_read_qc.csv", index=False)
            accepted = int(read_qc["accepted"].sum()) if len(read_qc) else 0
            rejected = int((~read_qc["accepted"]).sum()) if len(read_qc) else 0
            offsets = read_qc.loc[read_qc["accepted"], "offset"].value_counts().sort_index()
            orientations = read_qc.loc[read_qc["accepted"], "orientation"].value_counts().sort_index()
            qc_summary_rows.append({
                "locus": locus_name,
                "fragment": fragment,
                "condition": condition,
                "file": sample["fasta"],
                "raw_reads": len(records),
                "accepted_reads": accepted,
                "rejected_reads": rejected,
                "acceptance_rate": accepted / len(records) if records else 0.0,
                "mean_identity": read_qc.loc[read_qc["accepted"], "identity"].mean() if accepted else np.nan,
                "min_identity": read_qc.loc[read_qc["accepted"], "identity"].min() if accepted else np.nan,
                "max_identity": read_qc.loc[read_qc["accepted"], "identity"].max() if accepted else np.nan,
                "observed_offsets": "; ".join(f"{k}:{v}" for k, v in offsets.items()),
                "best_orientations": "; ".join(f"{k}:{v}" for k, v in orientations.items()),
                "ignored_non_CT_at_C_sites": int(sample_counts[condition]["ignored_non_CT_at_C_sites"].sum()),
            })
        if "WT" not in sample_counts or "Peloria" not in sample_counts:
            raise ValueError(f"Both WT and Peloria samples are required for {locus_name}/{fragment}")
        wt = sample_counts["WT"]
        pel = sample_counts["Peloria"]
        merged = ann.copy()
        merged["WT_methylated_reads"] = wt["methylated_reads"]
        merged["WT_total_reads"] = wt["total_reads"]
        merged["Peloria_methylated_reads"] = pel["methylated_reads"]
        merged["Peloria_total_reads"] = pel["total_reads"]
        merged = apply_min_coverage(merged, args.min_coverage)
        all_matrices.append(merged)

    matrix = pd.concat(all_matrices, ignore_index=True)
    matrix.to_csv(outdir / "cytosine_matrix.csv", index=False)
    pd.DataFrame(qc_summary_rows).to_csv(outdir / "alignment_qc.csv", index=False)
    metrics = add_epibit_metrics(matrix, alpha=args.alpha)
    summary = summarize_by_context(metrics)
    ranked, thresholds = dimensionality(metrics)
    metrics.to_csv(outdir / "derived_metrics.csv", index=False)
    summary.to_csv(outdir / "summary_by_context.csv", index=False)
    ranked.to_csv(outdir / "ranked_dimensionality.csv", index=False)
    thresholds.to_csv(outdir / "dimensionality_thresholds.csv", index=False)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="EpiBit cytosine-level epiallele analysis")
    sub = parser.add_subparsers(dest="mode", required=True)

    p1 = sub.add_parser("from-matrix", help="Compute metrics from a cytosine-level count matrix")
    p1.add_argument("--matrix", required=True, help="CSV matrix with WT/Peloria methylated and total read columns")
    p1.add_argument("--outdir", required=True)
    p1.add_argument("--alpha", type=float, default=ALPHA)
    p1.add_argument("--min-coverage", type=int, default=None, help="Optional coverage threshold in both conditions. Default: none.")
    p1.set_defaults(func=run_from_matrix)

    p2 = sub.add_parser("align", help="Align bisulfite reads and compute metrics")
    p2.add_argument("--config", required=True, help="JSON config listing loci, reference sequences and FASTA samples")
    p2.add_argument("--outdir", required=True)
    p2.add_argument("--min-identity", type=float, default=MIN_IDENTITY)
    p2.add_argument("--min-overlap", type=int, default=50)
    p2.add_argument("--alpha", type=float, default=ALPHA)
    p2.add_argument("--min-coverage", type=int, default=None, help="Optional post-alignment coverage threshold. Default: none.")
    p2.set_defaults(func=run_align)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
