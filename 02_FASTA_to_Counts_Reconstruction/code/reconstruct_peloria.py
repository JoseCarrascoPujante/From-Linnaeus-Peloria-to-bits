#!/usr/bin/env python3
from pathlib import Path
import argparse, csv
from epibit_core import load_fasta_dict, count_cytosines, context_of_cytosine, corrected_probability
from config import SAMPLES, ALPHA, MIN_IDENTITY, MIN_OVERLAP

def main():
    ap=argparse.ArgumentParser(description='Reconstruct Peloria/LvHIRZ cytosine count matrices from raw FASTA')
    ap.add_argument('--root', default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument('--outdir', default=None)
    args=ap.parse_args(); root=Path(args.root); out=Path(args.outdir or root/'outputs'/'peloria_reconstructed'); out.mkdir(parents=True,exist_ok=True)
    refs=load_fasta_dict(root/'reference_data'/'reference_sequences.fasta')
    by_locus={'Lcyc':[],'LvHIRZ':[]}; qc=[]
    for fragment,cfg in SAMPLES.items():
        ref=refs[fragment]
        condition_counts={}; condition_qc={}
        for cond in ('WT','Peloria'):
            counts,q=count_cytosines(root/'raw_data'/cfg[cond], ref, cfg['orientation'], MIN_IDENTITY, MIN_OVERLAP)
            condition_counts[cond]=counts; condition_qc[cond]=q
            qc.append({'locus':cfg['locus'],'fragment':fragment,'condition':cond,'file':cfg[cond],**q})
        for pos0,b in enumerate(ref):
            if b!='C': continue
            wt=condition_counts['WT'][pos0]; pe=condition_counts['Peloria'][pos0]
            if wt[1] < 1 or pe[1] < 1: continue
            by_locus[cfg['locus']].append({
                'fragment':fragment,'position':pos0+1,'context':context_of_cytosine(ref,pos0,fragment),
                'WT methylated':wt[0],'WT total':wt[1],
                'Peloria methylated':pe[0],'Peloria total':pe[1],
                'WT p':corrected_probability(wt[0],wt[1],ALPHA),
                'Peloria q':corrected_probability(pe[0],pe[1],ALPHA)})
    for locus,rows in by_locus.items():
        for i,r in enumerate(rows,1): r['site_id']=f'{locus}_C{i:03d}'
        cols=['site_id','fragment','position','context','WT methylated','WT total','Peloria methylated','Peloria total','WT p','Peloria q']
        with open(out/f'{locus}_cytosine_matrix_reconstructed.csv','w',newline='',encoding='utf-8') as fh:
            w=csv.DictWriter(fh,fieldnames=cols);w.writeheader();w.writerows(rows)
    with open(out/'alignment_QC_reconstructed.csv','w',newline='',encoding='utf-8') as fh:
        cols=['locus','fragment','condition','file','raw_reads','accepted_reads','rejected_reads','acceptance_rate','mean_identity','min_identity','max_identity','observed_offsets','orientations','ignored_non_CT_at_C_sites']
        w=csv.DictWriter(fh,fieldnames=cols);w.writeheader()
        for row in qc:
            row=dict(row)
            for k in ('observed_offsets','orientations'): row[k]='; '.join(f'{a}:{b}' for a,b in row[k].items())
            w.writerow(row)
    print(f'Wrote reconstructed matrices and QC to {out}')
if __name__=='__main__': main()
