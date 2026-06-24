#!/bin/usr python

import os
import pandas as pd
from Bio import SeqIO
import numpy as np
from itertools import groupby
import re
from operator import itemgetter
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import scipy as sp
from glob import glob
from collections import Counter
import matplotlib as mpl
from matplotlib.lines import Line2D
import random
from copy import deepcopy

# Third script to run for AAS Pipeline when using FragPipe for Val search
# Based on script from Tsour et al., Nature 2026 
# Last Updated 06-24-2026

"""
    This code reads in dictionaries generated in AA_subs_validation1.py
    Requires dataset metadata and result files from validation search.
    Output is MTP dict filtered for MTPs quantified in validation search 
    and for peptides with b/y ion evidence covering site of AAS.
"""

print("Beginning validation search 2 for FRAGPIPE data.")

### Set Directories
print("Setting directories...")
home_dir        = '/home/maropakis.a/'
scratch_dir     = '/scratch/maropakis.a/'
# code_dir      = home_dir + 'scripts/'

Frag_dir        = scratch_dir + 'Frag_outputs/results/        /  ' # edit to include data folder
aas_dir         = scratch_dir + 'AAS_Pipeline/  '  # edit to include data folder
database_dir    = scratch_dir + 'Dependencies/FASTA/'
sample_map = pd.read_excel(scratch_dir + 'Dependencies/sample_map/        .xlsx') # edit sample map

mtp_dict = pickle.load(open(aas_dir + 'qMTP_dict.p', 'rb'))
samples = list(mtp_dict.keys())

print("Directories loaded, beginning processing.")

""" Identify MTPs found in the validation FragPipe search """
print("Building validation evidence dict")
# initialize and save dictionary of validation evidence data 
val_evidence_dict = {} # create dict with all DP combined_ion.tsv data from all tmt sets 
for s in samples:
    print(s)
    evidence = pd.read_csv(Frag_dir+'*_1/psm.tsv', sep='\t', engine='python')
    # filter peptides for Probability (PEP) >= 0.99 and Purity (PIF) >= 0.8
    evidence = evidence.loc[(evidence['Probability']>=0.99) & (evidence['Purity']>=0.8),:]
    val_evidence_dict[s] = evidence
pickle.dump(val_evidence_dict, open(aas_dir+'Validation_search_evidence_dict.p', 'wb'))
#val_evidence_dict = pickle.load(open(aas_dir+'Validation_search_evidence_dict.p', 'rb'))

### determine if sequence is identified in validation search
def seq_in_val_search(idx, tmt_set):
    """
    Input: index of mtp_dict[tmt_set], tmt_set(sample)
    Output: if peptide is found, output = [index in mtp_dict, index in psm.tsv, index in mtp list at idx]
    """
    val_ev_df = val_evidence_dict[tmt_set]
    ev_seqs = list(val_ev_df['Peptide'].values)

    mtp_list = mtp_dict[tmt_set]['mistranslated sequence'][idx]
    # if >1 identified in val search, return as separate results
    for i, mtp in enumerate(mtp_list):
        if mtp in ev_seqs:
            return([idx, [i for i,x in enumerate(ev_seqs) if x==mtp], i])
        else:
            return None

### apply validation search function to each SAAP entry 
print("Identifying mtps that are found in regular search")
val_hit_lists = {s:[] for s in samples}
for s, s_dict in mtp_dict.items():
    for idx in s_dict['Raw file'].keys():
        result = seq_in_val_search(idx,s)
        if result:
            val_hit_lists[s].append(result)

### Loop through lists of results, create new dict of val SAAPs with link to psm.tsv file index
val_mtp_dict = {}
for s,val_list in val_hit_lists.items():
    print(s)
    val_mtp_dict[s] = {k:{} for k in mtp_dict[s].keys()}
    val_mtp_dict[s]['idx_val_evidence'] = {}

    for i, val in enumerate(val_list):
        mtp_idx = val[0]
        seq_idx = val[2]
        ev_idx = val[1]

        for k in mtp_dict[s].keys():
            if (isinstance(mtp_dict[s][k][mtp_idx], list)) and len(mtp_dict[s][k][mtp_idx])>0: # this is to make sure that we are extracting the correct AAS data and q-values for the mtp found out of list of mtps
                val_mtp_dict[s][k][i] = mtp_dict[s][k][mtp_idx][seq_idx]
            else:
                val_mtp_dict[s][k][i] = mtp_dict[s][k][mtp_idx]
        val_mtp_dict[s]['idx_val_evidence'][i] = ev_idx
    print(len(val_mtp_dict[s]['idx_val_evidence']))
pickle.dump(val_mtp_dict, open(aas_dir+'Validated_MTP_dict.p', 'wb'))
val_mtp_dict = pickle.load(open(aas_dir+'Validated_MTP_dict.p', 'rb'))

### determine the number of fragment ions supporting site of AAS
def n_frags_over_MTP(frag_match, mtp, sub_idx, tmt_set):
    """
    Input: fragment list from MSFragger tsv (nterm{n}=b{n}, cterm{n}=y{n}), peptide sequence, index of AAS on sequence, tmt_set
    Output: number of fragment ions covering site of AAS
    """
    count = 0
    L = len(mtp)
    for f, frag in enumerate(frag_match):
        mtp_frag=0
        ion = frag.split('^')[0] # drop charge state
        if 'nterm' in ion:
            frag_end = int(ion[5:])
            if frag_end>sub_idx:
                mtp_frag = 1
        elif 'cterm' in ion:
            frag_n = int(ion[5:])
            if L-frag_n <= sub_idx:
                mtp_frag=1
        count+=mtp_frag

    return(count)

### build dict of fragment tsvs keyed by raw file stem to avoid re-reading
frag_tsv_dict = {}
for f in glob(Frag_dir+'*_1/*FR*.tsv') + glob(Frag_dir+'*_1/*fraction*.tsv'):
    raw = os.path.splitext(os.path.basename(f))[0]
    frag_tsv_dict[raw] = pd.read_csv(f, sep='\t', low_memory=False)

### apply function and annotate val_mtp_dict
"""requires per-fraction MSFragger .tsv files (FragPipe output)"""
for s in samples:
    print(s)
    ev = val_evidence_dict[s]

    val_mtp_dict[s]['fragment_evidence'] = {}
    for k,v in val_mtp_dict[s]['aa subs'].items():
        seq = val_mtp_dict[s]['mistranslated sequence'][k]
        bp = val_mtp_dict[s]['DP Base Sequence'][k]
        sub_idx = [i for i,x in enumerate(bp) if seq[i]!=x][0]

        ev_idx = val_mtp_dict[s]['idx_val_evidence'][k]
        val_mtp_dict[s]['fragment_evidence'][k] = 0
        for idx in ev_idx:
            row = ev.iloc[idx,:]
            spectrum = row['Spectrum'] # FragPipe: RawStem.Scan.Scan.Charge
            raw_file = spectrum.rsplit('.',3)[0]
            scan = int(spectrum.split('.')[-3])
            frag_df = frag_tsv_dict.get(raw_file)
            if frag_df is None:
                continue
            scan_row = frag_df.loc[(frag_df['scannum']==scan) & (frag_df['hit_rank']==1),:]
            if len(scan_row)>0:
                matches_val = scan_row['fragments'].values[0]
                if isinstance(matches_val, str):
                    frag_match = matches_val.split(';')
                    frag_match = [x for x in frag_match if x] # drop trailing empty
                    count_frags = n_frags_over_MTP(frag_match, seq, sub_idx, s)
                    if count_frags>val_mtp_dict[s]['fragment_evidence'][k]:
                        val_mtp_dict[s]['fragment_evidence'][k] = count_frags

### filter val_mtp_dict for those with b/y ion evidence
val_ion_mtp_dict = {s:{} for s in samples}
for s in samples:
    ion_idx = [i for i,x in val_mtp_dict[s]['fragment_evidence'].items() if x>1]
    for k,v in val_mtp_dict[s].items():
        val_ion_mtp_dict[s][k] = {i:x for i,x in v.items() if i in ion_idx}
pickle.dump(val_ion_mtp_dict, open(aas_dir+'Ion_validated_MTP_dict.p', 'wb'))
print('Validation 2 completed successfully.')