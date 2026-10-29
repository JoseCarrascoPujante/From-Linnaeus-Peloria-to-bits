#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, io, zipfile
import pandas as pd

def read_member(zip_path, suffix):
    with zipfile.ZipFile(zip_path) as z:
        name=[n for n in z.namelist() if n.endswith(suffix)][0]
        return pd.read_csv(io.BytesIO(z.read(name)))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--root',default=str(Path(__file__).resolve().parents[1])); ap.add_argument('--outdir',default=None); args=ap.parse_args()
    root=Path(args.root); out=Path(args.outdir or root/'tables'); out.mkdir(parents=True,exist_ok=True)
    met=root/'outputs'/'epibit_metrics'
    lcyc=pd.read_csv(met/'Lcyc_context_summary.csv'); lv=pd.read_csv(met/'LvHIRZ_context_summary.csv')
    fwa=read_member(root.parent/'04_Supplementary_Data'/'Data_S3_FWA.zip','FWA_summary.csv')
    sup=read_member(root.parent/'04_Supplementary_Data'/'Data_S4_SUPERMAN.zip','SUPERMAN_summary.csv')
    h06=read_member(root.parent/'04_Supplementary_Data'/'Data_S6_H06.zip','H06_summary.csv')
    h06w=read_member(root.parent/'04_Supplementary_Data'/'Data_S6_H06.zip','H06_window_sensitivity.csv')
    arch=read_member(root.parent/'04_Supplementary_Data'/'Data_S7_architecture_space.zip','architecture_source_and_coordinates.csv')
    perm=read_member(root.parent/'04_Supplementary_Data'/'Data_S7_architecture_space.zip','context_label_permutation_summary.csv')

    # Main Table 1
    lt=lcyc[lcyc.context=='Total'].iloc[0]; lcg=lcyc.set_index('context'); dims=pd.read_csv(met/'Lcyc_dimensionality.csv')
    t1=pd.DataFrame([
        ['Empirical matrix','Informative cytosines',f"{int(lt.N)}: {int(lcg.loc['CG','N'])} CG, {int(lcg.loc['CHG','N'])} CHG, {int(lcg.loc['CHH','N'])} CHH"],
        ['Mean methylation','All Lcyc sites',f"WT {lt.WT_mean:.3f}; Peloria {lt.Peloria_mean:.3f}"],
        ['Directional divergence','KL(Peloria||WT)',f"{lt.KL_Peloria_WT:.2f} bits"],
        ['Symmetric separation','JSD / Dsym',f"{lt.JSD:.2f} / {lt.D_sym:.2f} bits"],
        ['Context allocation','CG / CHG / CHH',f"{lcg.loc['CG','KL_Peloria_WT']:.2f} / {lcg.loc['CHG','KL_Peloria_WT']:.2f} / {lcg.loc['CHH','KL_Peloria_WT']:.2f} bits"],
        ['CHG enrichment','E_CHG',f"{(lcg.loc['CHG','KL_Peloria_WT']/lt.KL_Peloria_WT)/(lcg.loc['CHG','N']/lt.N):.2f}x"],
        ['Non-CG dominance','CHG + CHH',f"{lcg.loc['CHG','KL_Peloria_WT']+lcg.loc['CHH','KL_Peloria_WT']:.2f} bits; {(lcg.loc['CHG','KL_Peloria_WT']+lcg.loc['CHH','KL_Peloria_WT'])/lt.KL_Peloria_WT*100:.1f}%"],
        ['Dimensionality','N50 / N80 / N90 / N95',' / '.join(str(int(x)) for x in dims.N_sites)+' sites'],
        ['Locus specificity','Lcyc vs LvHIRZ KL',f"{lt.KL_Peloria_WT:.2f} vs {lv[lv.context=='Total'].iloc[0].KL_Peloria_WT:.3f} bits"],
    ],columns=['Category','Metric','Result'])
    t1.to_csv(out/'Main_Table_1.csv',index=False)

    # Main Table 2 (Peloria pseudocounts) recomputed separately by make_figures; copy from archived original if present not needed here
    # Main Table 3
    pel=arch[arch.system=='Peloria'].iloc[0]
    rows=[]
    for sys,label in [('Peloria','Peloria / Lcyc'),('FWA','FWA'),('SUPERMAN','SUPERMAN / Wa-1'),('H06','H06 / IL8-3')]:
        r=arch[arch.system==sys].iloc[0]
        rows.append([label,int(r.N),r.KL_bits,r.KL_bits/r.N,r.q_CG,r.q_CHG,r.q_CHH,r.E_CG,r.E_CHG,r.E_CHH,int(r.N80)])
    pd.DataFrame(rows,columns=['System','N','KL_bits','KL_per_site','q_CG','q_CHG','q_CHH','E_CG','E_CHG','E_CHH','N80']).to_csv(out/'Main_Table_3.csv',index=False)
    arch[['system','X_nonCG_dominance','Y_log2_ECHG','context_breadth_B','R80','Z_log2_ECG','ILR1_CG_vs_nonCG','ILR2_CHG_vs_CHH']].to_csv(out/'Main_Table_4_architecture_coordinates.csv',index=False)

    # Supplementary exact tables
    fwa.to_csv(out/'Supplementary_Table_S10_FWA_exact.csv',index=False)
    sup.to_csv(out/'Supplementary_Table_S11_SUPERMAN_exact.csv',index=False)
    h06w.to_csv(out/'Supplementary_Table_S12_H06_window_exact.csv',index=False)
    arch.to_csv(out/'Supplementary_Table_S13_architecture_exact.csv',index=False)
    perm.to_csv(out/'Supplementary_context_permutation_exact.csv',index=False)
    print(f'Wrote manuscript/source tables to {out}')
if __name__=='__main__': main()
