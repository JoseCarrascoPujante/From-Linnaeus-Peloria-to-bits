#!/usr/bin/env python3
"""Cross-check of epibit.py against an independent implementation.

Both start from the reconstructed Lcyc/LvHIRZ count matrices; agreement of regional
totals, per-context directional KL and dimensionality thresholds shows that
both compute the same EpiBit quantities.
"""
from pathlib import Path
import argparse, sys
import numpy as np, pandas as pd

ap = argparse.ArgumentParser(); ap.add_argument('--root', required=True); r = Path(ap.parse_args().root)
sys.path.insert(0, str(r/'02_FASTA_to_Counts_Reconstruction'/'code'))
from epibit_core import corrected_probability, site_metrics  # noqa: E402
from config import ALPHA  # noqa: E402

TOL = 1e-9
ok = True
def check(label, a, b, tol=TOL):
    global ok
    good = (a == b) if tol == 0 else abs(float(a) - float(b)) <= tol
    print('PASS' if good else 'FAIL', label, a, b); ok &= good

for locus in ['Lcyc', 'LvHIRZ']:
    m = pd.read_csv(r/'02_FASTA_to_Counts_Reconstruction'/'outputs'/'peloria_reconstructed'/f'{locus}_cytosine_matrix_reconstructed.csv')
    p = [corrected_probability(a, b, ALPHA) for a, b in zip(m['WT methylated'], m['WT total'])]
    q = [corrected_probability(a, b, ALPHA) for a, b in zip(m['Peloria methylated'], m['Peloria total'])]
    s = pd.DataFrame([site_metrics(x, y) for x, y in zip(p, q)]); s['context'] = m.context
    h = r/'05_Reproduction_Outputs'/'historical_epibit'/locus
    region = pd.read_csv(h/'region_summary.csv').iloc[0]
    check(f'{locus} N', int(region.N_cytosines), len(m), 0)
    check(f'{locus} WT_mean', region.WT_mean_p, np.mean(p)); check(f'{locus} Peloria_mean', region.Peloria_mean_q, np.mean(q))
    check(f'{locus} KL', region.KL_Peloria_WT_bits, s.KL_Peloria_WT.sum())
    check(f'{locus} Dsym', region.Dsym_bits, s.D_sym.sum()); check(f'{locus} JSD', region.JSD_bits, s.JSD.sum())
    ctx = pd.read_csv(h/'context_summary.csv').set_index('context')
    for c in ctx.index:
        check(f'{locus} {c} KL', ctx.loc[c, 'KL_Peloria_WT_bits'], s.KL_Peloria_WT[s.context == c].sum())
    cum = np.cumsum(np.sort(s.D_sym.to_numpy())[::-1])/s.D_sym.sum()
    ind = [int(np.argmax(cum >= t)) + 1 for t in (0.5, 0.8, 0.9, 0.95)]
    check(f'{locus} dimensionality', [int(x) for x in pd.read_csv(h/'dimensionality_thresholds.csv').k_sites], ind, 0)
if not ok:
    sys.exit(1)
print('HISTORICAL EPIBIT VALIDATION PASSED')
