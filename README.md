# From Linnaeus' Peloria to bits

Code and data for the manuscript "From Linnaeus' Peloria to bits: a CHG-dominated architecture of natural epiallelic divergence".

## Contents

| Folder | Contents |
|---|---|
| `01_EpiBit_Code/` | EpiBit analysis program (`code/epibit.py`), its unit test and its input template (`data/cytosine_matrix_template.csv`) |
| `02_FASTA_to_Counts_Reconstruction/` | Targeted bisulfite FASTA reads (`raw_data/`), reference fragments (`reference_data/`) and the code that converts them into the Lcyc and LvHIRZ cytosine count matrices (`outputs/`) |
| `03_Comparative_Extension/` | Code that produces the manuscript tables and data-driven figures (Peloria, FWA, SUPERMAN, H06 and architecture space), with its outputs |
| `04_Supplementary_Data/` | Supplementary Data S1–S7 |
| `05_Reproduction_Outputs/` | Results of `epibit.py` on the reconstructed matrices |
| `06_EpiBit_Benchmark/` | Benchmark of EpiBit against global and context-specific mean methylation change (Figure 7) |
| `workflow/` | `run_all.py` and the scripts it calls |

## Requirements

Python 3.10 or later with numpy, pandas, openpyxl and matplotlib (tested with Python 3.10.12, numpy 2.2.6, pandas 2.3.3, openpyxl 3.1.5 and matplotlib 3.10.9).

## Reproduction

```bash
python workflow/run_all.py
```

The workflow runs the unit test of `epibit.py`; converts the FASTA reads into count matrices; checks that Supplementary Data S1, S2 and the Peloria rows of S7 match those matrices (`python workflow/build_data_packages.py --write` rebuilds them); runs `epibit.py` at α = 0.5 and its pseudocount sensitivity analysis; cross-checks `epibit.py` against an independent implementation; and produces the comparative tables and figures. A successful run ends with `ALL END-TO-END REPRODUCTION STEPS COMPLETED SUCCESSFULLY`. The benchmark runs separately:

```bash
python 06_EpiBit_Benchmark/benchmark_epibit_vs_simple.py
```

Expected results: Lcyc, 98 informative cytosines, directional KL 87.418 bits, Dsym 83.081 bits, JSD 16.315 bits, N50/N80/N90/N95 = 12/25/34/43; LvHIRZ, 29 informative cytosines, directional KL 0.934 bits.

## Provenance

`epibit.py` is the EpiBit program of the original submission. The FASTA-to-count step follows the manuscript methods: each read is evaluated at every ungapped offset with at least 50 aligned bases, in the read orientation of its fragment (both orientations for LvHIRZ), and accepted when its best bisulfite-aware identity is at least 0.85. Compared with v1.0.0 of this repository (https://doi.org/10.5281/zenodo.22260591), the LvHIRZ matrix is identical and two additional Lcyc fragment-3 reads are accepted, which adds one read to 20 of the 98 cytosines. All other code is new.

## Data, license and citation

The FASTA reads were provided by A. E. Ehrenhofer-Murray (Ehrenhofer-Murray, 2021); accessions of the public datasets are given in the manuscript. Code is released under the MIT License and data under CC-BY 4.0, with the third-party terms stated in `LICENSE`. Citation metadata are in `CITATION.cff`.
