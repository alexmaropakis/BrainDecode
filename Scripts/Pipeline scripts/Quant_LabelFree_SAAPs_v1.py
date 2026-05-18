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

# Fourth and final script to run for label free AAS Pipeline
# Last updated 05-07-2026 by Alex Maropakis


""" This script is for use with label-free datasets. It follows the 
    AAS_detection and AAS_validation scripts. 
    It takes in highly confidently identified peptides with AAS, 
    and quantifies them. Output metrics include normalized abundances, 
    and ratios of peptide with AAS to its base peptide, at precursor ion levels, 
    stored in a dict that structured by unique MTP-BP pairs.
    This script also generates a dictionary of quant data organized by AAS type.
"""

print("Beginning final quantification...")

### Set Directories
print("Setting directories...")
home_dir            = '/home/maropakis.a/'
scratch_dir         = '/scratch/maropakis.a/'
# code_dir          = home_dir + 'scripts/'

MQ_dir              = scratch_dir + 'MQ_outputs/    /Val/combined/txt/' # edit to include data folder
aas_dir             = scratch_dir + 'AAS_Pipeline/  '                  # edit to include data folder

mtp_dict            = pickle.load(open(aas_dir+'Ion_validated_MTP_dict.p', 'rb'))
samples             = list(mtp_dict.keys())
val_evidence_dict   = pickle.load(open(aas_dir+'Validation_search_evidence_dict.p', 'rb'))

print("Initializing...")
unq_pairs = []
unq_pair_dicts = []
for s, s_dict in mtp_dict.items():
    print(s)
    for i, v in s_dict['aa subs'].items():
        sub = v
        mtp = s_dict['mistranslated sequence'][i]
        bp  = s_dict['DP Base Sequence'][i]
        ev_idx = s_dict['idx_val_evidence'][i]
        pair = [mtp, bp]
        if pair not in unq_pairs:
            unq_pairs.append(pair)
            unq_pair_dicts.append({'MTP':mtp, 'BP':bp, 'AAS':sub,
                                'tissue_list':[s], 'ev_idx_list':ev_idx})
            unq_pair_dicts[len(unq_pair_dicts)-1]['tissue_evidence'] = {s:[ev_idx]}
        else:
            curr_dict_idx = [j for j, x in enumerate(unq_pair_dicts)
                            if (x['MTP']==mtp) and (x['BP']==bp)][0]
            unq_pair_dicts[curr_dict_idx]['tissue_list'].append(s)
            [unq_pair_dicts[curr_dict_idx]['ev_idx_list'].append(x) for x in ev_idx]
            if s in unq_pair_dicts[curr_dict_idx]['tissue_evidence'].keys():
                [unq_pair_dicts[curr_dict_idx]['tissue_evidence'][s].append(x) for x in ev_idx]
            else:
                unq_pair_dicts[curr_dict_idx]['tissue_evidence'][s] = ev_idx
print(f"Unique MTP-BP pairs:  {len(unq_pairs)}")

MTP_quant_dict = {}
for i, pair in enumerate(unq_pair_dicts):

    curr_dict = {'MTP_seq':pair['MTP'], 
        'BP_seq':pair['BP'], 
        'sub_index':[i for i,x in enumerate(pair['MTP']) if pair['MTP'][i]!=pair['BP'][i]], 
        'aa_sub':pair['AAS'],'tissues':list(set(pair['tissue_list'])), 
        'tissue_evidence':pair['tissue_evidence'],
        'MTP_PrecInt':{x:np.nan for x in samples}, 
        'BP_PrecInt':{x:np.nan for x in samples},
        'Prec_ratio':{x:np.nan for x in samples}, 
        'Norm_MTP_PrecInt':{x:np.nan for x in samples}, 
        'Norm_BP_PrecInt':{x:np.nan for x in samples}
    }
    MTP_quant_dict[i] = curr_dict

def median_normalize(sample_raw):
    """
    Input: list of raw precursor intensities for tissue
    Output: median-normalized list of precursor intensities for tissue
    """
    sample_median = np.median(sample_raw)
    sample_norm = [x/sample_median for x in sample_raw]
    return(sample_norm)

def precursor_intensity_quant(k, tissue):
    """
    Input: k=index of MTP in MTP_quant_dict, tissue
    Output: precursor intensity of mtp, its bp, their ratio
            Precursor intensities are log-transformed values of those reported in evidence.txt result file from validation MQ search.
            Precursor intensities represent the sum of all precursors mapped to the peptide sequence (across multiple fractions, charge states)
    """
    bp = MTP_quant_dict[k]['BP_seq']
    mtp = MTP_quant_dict[k]['MTP_seq']
    ev_df = val_evidence_dict[tissue]
    bp_ev_df = ev_df.loc[ev_df['Sequence']==bp,:]
    bp_prec_int = np.sum([x for x in bp_ev_df['Intensity'].values if ~np.isnan(x)])
    norm_bp_prec_int = np.sum([x for x in bp_ev_df['Intensity'].values if ~np.isnan(x)])
    
    mtp_ev_df = ev_df.loc[ev_df['Sequence']==mtp,:]
    mtp_prec_int = np.sum([x for x in mtp_ev_df['Intensity'].values if ~np.isnan(x)])
    norm_mtp_prec_int = np.sum([x for x in mtp_ev_df['Intensity'].values if ~np.isnan(x)])
  
    prec_ratio = np.log2(mtp_prec_int/bp_prec_int)
    return([mtp_prec_int, bp_prec_int, norm_mtp_prec_int, norm_bp_prec_int, prec_ratio])




print("Quantifying SAAPS")
for s in samples:
    print(s)
    tissue=s
    ev_df = val_evidence_dict[s]
    for col in [x for x in ev_df if x=='Intensity']:
        ev_df[col+'_Normalized'] = median_normalize(ev_df[col].values)
    
    for k,v in MTP_quant_dict.items():
        if s in v['tissues']:
            mtp = v['MTP_seq']
            bp = v['BP_seq']
            bp_ev = ev_df.loc[ev_df['Sequence']==bp,:]
            mtp_ev = ev_df.loc[ev_df['Sequence']==mtp,:]

            mtp_prec_int, bp_prec_int, norm_mtp_prec_int, norm_bp_prec_int, prec_ratio = precursor_intensity_quant(k, tissue)
            v['BP_PrecInt'][tissue] = bp_prec_int
            v['MTP_PrecInt'][tissue] = mtp_prec_int
            v['Norm_BP_PrecInt'][tissue] = norm_bp_prec_int
            v['Norm_MTP_PrecInt'][tissue] = norm_mtp_prec_int
            v['Prec_ratio'][tissue] = prec_ratio

print("Saving results to MTP_quant_dict.p")
pickle.dump(MTP_quant_dict, open(aas_dir+'MTP_quant_dict.p', 'wb'))

print("Quantification of label free SAAPs finished!")