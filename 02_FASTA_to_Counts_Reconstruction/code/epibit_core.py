#!/usr/bin/env python3
"""Core EpiBit routines used by the reproducible Peloria workflow.

FASTA-to-count step implemented from the manuscript methods: ungapped
bisulfite-aware alignment of each raw read to its reference fragment, read
acceptance at identity >= 0.85 over >= 50 aligned bases, and cytosine-level
methylated/informative counts from the accepted reads.
"""
from __future__ import annotations
import math
from pathlib import Path
from collections import defaultdict, Counter

import numpy as np

COMP = str.maketrans("ACGTN", "TGCAN")

def reverse_complement(seq: str) -> str:
    return seq.upper().translate(COMP)[::-1]

def read_fasta(path):
    name = None; seq=[]
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line=line.strip()
            if not line: continue
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(seq).upper()
                name=line[1:]; seq=[]
            else:
                seq.append(line)
    if name is not None:
        yield name, ''.join(seq).upper()

def load_fasta_dict(path):
    return {name.split()[0]:seq for name,seq in read_fasta(path)}

def context_of_cytosine(ref: str, pos0: int, fragment: str | None = None) -> str:
    """Classify a forward-reference C as CG/CHG/CHH.

    The two terminal C positions in LCYC_peloria_fragment4 lack enough
    downstream sequence in the deposited fragment to infer context de novo.
    They are retained as CHH to match the deposited annotation.
    """
    if ref[pos0] != 'C':
        raise ValueError('position is not C')
    if pos0 + 1 >= len(ref):
        return 'CHH' if fragment == 'LCYC_peloria_fragment4' else 'C_unknown'
    b1=ref[pos0+1]
    if b1 == 'G': return 'CG'
    if pos0 + 2 >= len(ref):
        return 'CHH' if fragment == 'LCYC_peloria_fragment4' else 'C_unknown'
    b2=ref[pos0+2]
    if b1 in 'ACT' and b2 == 'G': return 'CHG'
    if b1 in 'ACT' and b2 in 'ACT': return 'CHH'
    return 'C_unknown'

def _one_hot(seq: str):
    codes = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
    return {b: (codes == ord(b)).astype(np.float64) for b in 'ACGT'}

def bisulfite_scan(seq: str, ref_onehot, ref_len: int, min_overlap: int):
    """Bisulfite-compatible matches of seq at every ungapped offset on the reference.

    Offset d places read position 0 at reference position d (negative offsets start
    the read upstream of the fragment). A reference C matches a read C (methylated)
    or T (converted); all other bases require identity. Only offsets with at least
    min_overlap aligned bases are returned: the acceptance rule itself bounds the
    search, so no offset has to be specified in advance.
    """
    q = _one_hot(seq); n = len(seq)
    matches = sum(np.correlate(ref_onehot[b], q[b], mode='full') for b in 'ACGT')
    matches = matches + np.correlate(ref_onehot['C'], q['T'], mode='full')
    offsets = np.arange(-(n - 1), ref_len)
    overlap = np.minimum(n - np.maximum(0, -offsets), ref_len - np.maximum(0, offsets))
    keep = overlap >= min_overlap
    return offsets[keep], np.rint(matches[keep]).astype(np.int64), overlap[keep]

def best_alignment(read, ref_onehot, ref_len, orientations, min_overlap):
    """Highest-identity alignment over all offsets and allowed orientations.

    Ties are broken by longer overlap, then smaller |offset|, then forward orientation.
    """
    best = None
    for ori in orientations:
        seq = read if ori == 'forward' else reverse_complement(read)
        off, m, ov = bisulfite_scan(seq, ref_onehot, ref_len, min_overlap)
        if off.size == 0:
            continue
        ident = m / ov
        cand = np.flatnonzero(ident == ident.max())
        cand = cand[ov[cand] == ov[cand].max()]
        cand = cand[np.abs(off[cand]) == np.abs(off[cand]).min()]
        i = cand[0]
        key = (float(ident[i]), int(ov[i]), -abs(int(off[i])), ori == 'forward')
        if best is None or key > best[0]:
            best = (key, int(off[i]), seq, ori)
    return best

def count_cytosines(fasta_path, ref, orientation, min_identity=0.85, min_overlap=50):
    """Methylated/informative counts per reference cytosine and read-level QC.

    Every QC statistic, including the offset histogram, refers to accepted reads.
    """
    orientations = orientation if isinstance(orientation, (list, tuple)) else [orientation]
    ref_onehot = _one_hot(ref)
    counts = defaultdict(lambda: [0, 0])  # methylated, total informative
    raw = accepted = ignored = 0
    identities = []; offset_counts = Counter(); orientation_counts = Counter()
    for _, read in read_fasta(fasta_path):
        raw += 1
        best = best_alignment(read, ref_onehot, len(ref), orientations, min_overlap)
        if best is None or best[0][0] < min_identity:
            continue
        (ident, ov, _, _), off, seq, ori = best
        accepted += 1; identities.append(ident)
        offset_counts[off] += 1; orientation_counts[ori] += 1
        rs = max(0, off); qs = max(0, -off)
        for j in range(ov):
            if ref[rs + j] != 'C':
                continue
            b = seq[qs + j]
            if b in ('C', 'T'):
                counts[rs + j][1] += 1
                if b == 'C':
                    counts[rs + j][0] += 1
            else:
                # Non-C/T calls (including N) are neither methylated nor informative.
                ignored += 1
    return counts, {
        'raw_reads': raw,
        'accepted_reads': accepted,
        'rejected_reads': raw - accepted,
        'acceptance_rate': accepted / raw if raw else 0.0,
        'mean_identity': sum(identities) / len(identities) if identities else float('nan'),
        'min_identity': min(identities) if identities else float('nan'),
        'max_identity': max(identities) if identities else float('nan'),
        'observed_offsets': dict(sorted(offset_counts.items())),
        'orientations': dict(sorted(orientation_counts.items())),
        'ignored_non_CT_at_C_sites': ignored,
    }

def corrected_probability(methylated, total, alpha=0.5):
    return (methylated + alpha)/(total + 2*alpha)

def entropy(p):
    if p <= 0 or p >= 1: return 0.0
    return -p*math.log2(p)-(1-p)*math.log2(1-p)

def kl_bernoulli(q,p):
    return q*math.log2(q/p)+(1-q)*math.log2((1-q)/(1-p))

def jsd_bernoulli(p,q):
    m=(p+q)/2
    return 0.5*kl_bernoulli(p,m)+0.5*kl_bernoulli(q,m)

def site_metrics(p,q):
    kqp=kl_bernoulli(q,p); kpq=kl_bernoulli(p,q)
    return {
        'H_WT':entropy(p),'H_Peloria':entropy(q),
        'KL_Peloria_WT':kqp,'KL_WT_Peloria':kpq,
        'D_sym':0.5*(kqp+kpq),'JSD':jsd_bernoulli(p,q)
    }
