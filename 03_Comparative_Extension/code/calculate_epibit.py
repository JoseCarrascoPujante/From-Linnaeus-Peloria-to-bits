#!/usr/bin/env python3
from pathlib import Path
import argparse, csv
from collections import defaultdict
from epibit_core import site_metrics

def read_csv(path):
    with open(path,encoding='utf-8') as fh: return list(csv.DictReader(fh))
def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1])); ap.add_argument('--indir',default=None); ap.add_argument('--outdir',default=None); args=ap.parse_args()
    root=Path(args.root); indir=Path(args.indir or root/'outputs'/'peloria_reconstructed'); out=Path(args.outdir or root/'outputs'/'epibit_metrics'); out.mkdir(parents=True,exist_ok=True)
    for locus in ['Lcyc','LvHIRZ']:
        rows=read_csv(indir/f'{locus}_cytosine_matrix_reconstructed.csv'); site=[]
        for r in rows:
            p=float(r['WT p']); q=float(r['Peloria q']); m=site_metrics(p,q); site.append({**r,**m})
        cols=list(site[0].keys())
        with open(out/f'{locus}_site_metrics.csv','w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=cols);w.writeheader();w.writerows(site)
        contexts=defaultdict(list)
        for r in site: contexts[r['context']].append(r)
        summary=[]
        order=[c for c in ['CG','CHG','CHH','C_unknown'] if c in contexts]
        for c in order:
            rr=contexts[c]; summary.append({'context':c,'N':len(rr),'WT_mean':sum(float(x['WT p']) for x in rr)/len(rr),'Peloria_mean':sum(float(x['Peloria q']) for x in rr)/len(rr),**{k:sum(float(x[k]) for x in rr) for k in ['H_WT','H_Peloria','KL_Peloria_WT','KL_WT_Peloria','D_sym','JSD']}})
        rr=site; summary.append({'context':'Total','N':len(rr),'WT_mean':sum(float(x['WT p']) for x in rr)/len(rr),'Peloria_mean':sum(float(x['Peloria q']) for x in rr)/len(rr),**{k:sum(float(x[k]) for x in rr) for k in ['H_WT','H_Peloria','KL_Peloria_WT','KL_WT_Peloria','D_sym','JSD']}})
        with open(out/f'{locus}_context_summary.csv','w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=list(summary[0].keys()));w.writeheader();w.writerows(summary)
        ranked=sorted(site,key=lambda x:float(x['D_sym']),reverse=True); tot=sum(float(x['D_sym']) for x in ranked); cum=0.0; thresholds={}
        for i,r in enumerate(ranked,1):
            cum += float(r['D_sym'])
            for t in (0.5,0.8,0.9,0.95):
                if t not in thresholds and cum/tot >= t: thresholds[t]=i
        with open(out/f'{locus}_dimensionality.csv','w',newline='',encoding='utf-8') as fh:
            w=csv.writer(fh);w.writerow(['threshold','N_sites','fraction_of_sites'])
            for t in (0.5,0.8,0.9,0.95):w.writerow([t,thresholds[t],thresholds[t]/len(ranked)])
    print(f'Wrote EpiBit metrics to {out}')
if __name__=='__main__':main()
