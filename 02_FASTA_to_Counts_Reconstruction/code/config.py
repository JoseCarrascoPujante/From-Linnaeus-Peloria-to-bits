"""Analysis configuration for the Peloria targeted amplicons.

Read orientation per fragment follows the manuscript methods (reads evaluated in
the appropriate orientation); alignment offsets are not configured because every
ungapped offset with at least MIN_OVERLAP aligned bases is evaluated.
"""
ALPHA = 0.5
MIN_IDENTITY = 0.85
MIN_OVERLAP = 50
SEED = 20261022
PERMUTATIONS = 200000

SAMPLES = {
    'LCYC_F3_R3': {
        'locus':'Lcyc', 'orientation':'forward',
        'WT':'LCYC_FR_8.fasta', 'Peloria':'LCYC_FR_9.fasta'},
    'LCYC_peloria_fragment3': {
        'locus':'Lcyc', 'orientation':'forward',
        'WT':'LCYC_fr3_8.fasta', 'Peloria':'LCYC_fr3_9.fasta'},
    'LCYC_peloria_fragment4': {
        'locus':'Lcyc', 'orientation':'reverse',
        'WT':'LCYC_fr4_8.fasta', 'Peloria':'LCYC_fr4_9.fasta'},
    'LvHIRZ_2': {
        'locus':'LvHIRZ', 'orientation':'reverse',
        'WT':'LvHIRZ_8.fasta', 'Peloria':'LvHIRZ_9.fasta'},
}
