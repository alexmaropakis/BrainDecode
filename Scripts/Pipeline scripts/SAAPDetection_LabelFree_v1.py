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
import multiprocessing as mp

# First script in AAS Pipeline for LABEL FREE datasets
#  Last updated 05-07-2026 by Alex Maropakis

### Set Directories
print("Setting directories...")
home_dir            = '/home/maropakis.a/'
scratch_dir         = '/scratch/maropakis.a/'
# code_dir          = home_dir + 'scripts/'
MQ_dir              = scratch_dir + 'MQ_outputs/    /DP/combined/txt/'          # edit to include data folder
aas_dir             = scratch_dir + 'AAS_Pipeline/  '                           # edit to include data folder
fasta_path          = scratch_dir + 'Dependencies/FASTA/HUMAN_GENOME.fna'       # edit for species
frameshift_dir      = scratch_dir + 'Dependencies/frame_translations/human/'    # edit for species
mod_df              = pd.read_excel('/scratch/maropakis.a/AAS_Pipeline/Modifications_table_091520.xlsx')
os.makedirs(aas_dir, exist_ok=True)
os.makedirs(frameshift_dir, exist_ok=True)

# read in list of tissues/experiments/samples for which you have independent MQ results
tissues             = open(scratch_dir+'Dependencies/sample_map/    .txt', 'r').read().split('\n')

### Amino acid molecular weight dictionary 
print(f"Assembling AAS mass lookup table...")
MW_dict = {
    "G":57.02147,"A":71.03712,"S":87.03203,"P":97.05277,"V":99.06842,
    "T":101.04768,"I":113.08407,"L":113.08407,"N":114.04293,"D":115.02695,
    "Q":128.05858,"K":357.257902,"E":129.0426,"M":131.04049,"H":137.05891,
    "F":147.06842,"R":156.10112,"C":160.030654,"Y":163.0633,"W":186.07932,
}
AAs = 'ACDEFGHIKLMNPQRSTVWY' # list of amino acids 

raw_subs = {f'{i} to {j}': MW_dict[j]-MW_dict[i] for i in MW_dict for j in MW_dict if i!=j}
raw_subs.pop('L to I', None); raw_subs.pop('I to L', None) # update dict to account for I and L having same mass
merged = {}
for k, v in raw_subs.items(): # unifies I and L
    merged[k.replace('I','L') if 'I' in k else k] = v
# Group by origin AA; store as numpy arrays for vectorized mass check
subs_dict      = {}   # {aa: [sub_strings]}
subs_mass_dict = {}   # {aa: np.array of mass shifts}
for a in AAs:
    items = [(k, v) for k, v in merged.items() if k.startswith(a)]
    subs_dict[a]      = [k for k,_ in items]
    subs_mass_dict[a] = np.array([v for _,v in items], dtype=np.float64)

def codonify(seq):
    """
    input: a nucleotide sequence (not necessarily a string)
    output: a list of codons
    """
    seq = str(seq)
    if len(seq) % 3:
        seq += 'N' * (3 - len(seq) % 3)
    return [seq[i:i+3] for i in range(0, len(seq), 3)]

def suffix_array(text, _step=16):
    """Analyze all common strings in the text.
    
    Short substrings of the length _step a are first pre-sorted. Then are the
    results repeatedly merged so that the guaranteed number of compared
    characters bytes is doubled in every iteration until all substrings are
    sorted exactly.
    
    Arguments:
        text:  The text to be analyzed.
        _step: Is only for optimization and testing. It is the optimal length
               of substrings used for initial pre-sorting. The bigger value is
               faster if there is enough memory. Memory requirements are
               approximately (estimate for 32 bit Python 3.3):
                   len(text) * (29 + (_size + 20 if _size > 2 else 0)) + 1MB
    
    Return value:      (tuple)
      (sa, rsa, lcp)
        sa:  Suffix array                  for i in range(1, size):
               assert text[sa[i-1]:] < text[sa[i]:]
        rsa: Reverse suffix array          for i in range(size):
               assert rsa[sa[i]] == i
        lcp: Longest common prefix         for i in range(1, size):
               assert text[sa[i-1]:sa[i-1]+lcp[i]] == text[sa[i]:sa[i]+lcp[i]]
               if sa[i-1] + lcp[i] < len(text):
                   assert text[sa[i-1] + lcp[i]] < text[sa[i] + lcp[i]]
    suffix_array(text='banana')
    ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [0, 1, 3, 0, 0, 2])
    
    Explanation: 'a' < 'ana' < 'anana' < 'banana' < 'na' < 'nana'
    The Longest Common String is 'ana': lcp[2] == 3 == len('ana')
    It is between  tx[sa[1]:] == 'ana' < 'anana' == tx[sa[2]:]
    """
    size  = len(text)
    step  = min(max(_step, 1), size)
    sa    = list(range(size))
    sa.sort(key=lambda i: text[i:i+step])
    rsa   = [0]*size
    grpstart = [False]*size + [True]
    igrp, stgrp = 0, None
    for i, pos in enumerate(sa):
        st = text[pos:pos+step]
        if st != stgrp:
            grpstart[igrp] = (igrp < i-1); stgrp = st; igrp = i
        rsa[pos] = igrp; sa[i] = pos
    grpstart[igrp] = (igrp < size-1)
    while grpstart.index(True) < size:
        nextgr = grpstart.index(True)
        while nextgr < size:
            igrp = nextgr
            nextgr = grpstart.index(True, igrp+1)
            glist = []
            for ig in range(igrp, nextgr):
                pos = sa[ig]
                if rsa[pos] != igrp: break
                glist.append((rsa[pos+step] if pos+step < size else -1, pos))
            glist.sort()
            for _, g in groupby(glist, key=itemgetter(0)):
                g = [x[1] for x in g]
                sa[igrp:igrp+len(g)] = g
                grpstart[igrp] = (len(g) > 1)
                for pos in g: rsa[pos] = igrp
                igrp += len(g)
        step *= 2
    return sa

### Build or load frame translations
def build_ref_translation(fasta_path, frame):
    """ build reference fasta for removal of homologous sequences """
    print(f"Generating translation for frame {frame}...")
    seqs = []
    for rec in SeqIO.parse(fasta_path, 'fasta'):
        seq = rec.seq.upper()
        fs = seq[frame-1:] if frame <= 3 else seq.reverse_complement()[(frame-3)-1:]
        seqs.append(str(seq(''.join(codonify(fs))).translate(to_stop=True)))
    W   = ''.join(seqs).replace('I', 'L')
    sa  = np.array(suffix_array(W), dtype=np.int32)
    pickle.dump(W,  open(os.path.join(frameshift_dir, f'W{frame}_aa_ambig.p'), 'wb'))
    pickle.dump(sa, open(os.path.join(frameshift_dir, f's{frame}a_ambig.p'), 'wb'))
    return W, sa

def load_frame(frame):
    wp = os.path.join(frameshift_dir, f'W{frame}_aa_ambig.p')
    sp = os.path.join(frameshift_dir, f's{frame}a_ambig.p')
    if os.path.exists(wp) and os.path.exists(sp):
        print(f"Loading frameshift translation for frame {frame}...")
        W  = pickle.load(open(wp, 'rb'))
        sa = pickle.load(open(sp, 'rb'))
        # Ensure numpy int32 (old files may be plain lists)
        if not isinstance(sa, np.ndarray):
            sa = np.array(sa, dtype=np.int32)
            pickle.dump(sa, open(sp, 'wb'))
        return W, sa
    return build_ref_translation(fasta_path, frame)

# Combine all 6 frame translations into one string separated by '|' to perform 3 SA searches per peptide instead of 18 
FRAME_SEP = '|'
combined_W_path  = os.path.join(frameshift_dir, 'W_combined_ambig.p')
combined_sa_path = os.path.join(frameshift_dir, 'sa_combined_ambig.p')

frames = {f: load_frame(f) for f in range(1, 7)}

if os.path.exists(combined_W_path) and os.path.exists(combined_sa_path):
    print("Loading combined genome string from disk...")
    W_combined  = pickle.load(open(combined_W_path, 'rb'))
    sa_combined = pickle.load(open(combined_sa_path, 'rb'))
    if not isinstance(sa_combined, np.ndarray):
        sa_combined = np.array(sa_combined, dtype=np.int32)
else:
    print(f"Building combined genome string...")
    W_combined  = FRAME_SEP.join(frames[f][0] for f in range(1, 7))
    sa_combined = np.array(suffix_array(W_combined), dtype=np.int32)
    pickle.dump(W_combined,  open(combined_W_path,  'wb'))
    pickle.dump(sa_combined, open(combined_sa_path, 'wb'))
    print(f"Combined genome string saved.")
# Keep individual frames for prot_nterm / prot_cterm
W_aa_ambiguous, sa_ambiguous = frames[1]

def SA_search(P, W, sa): 
    """ search for peptide in genome in suffix array
    Input: P = peptide sequence, W = translated genome, sa = sorted substrings of proteome
    Output: list of indices of string match
    """
    lp = len(P)
    n  = len(sa)
    lo, hi = 0, n
    # Left bound
    while lo < hi:
        mid = (lo + hi) >> 1
        a   = int(sa[mid])
        if W[a:a+lp] < P: lo = mid + 1
        else:              hi = mid
    s  = lo
    hi = n
    # Right bound
    while lo < hi:
        mid = (lo + hi) >> 1
        a   = int(sa[mid])
        if W[a:a+lp] <= P: lo = mid + 1
        else:               hi = mid
    return sa[s:lo]

# find homologous peptides
def _homology_worker(P): 
    """
    Gets a peptide and returns whether it has homolegous translation in the genome.
    """
    return P, any(
        len(SA_search(pfx + P, W_combined, sa_combined)) > 0
        for pfx in ('K', 'R', '*')
    )

def batch_homology_search(seq_list):
    """
    Search homology for a deduplicated list of sequences.
    Runs sequentially since this is called from within a pool worker.
    Returns dict {seq: bool}.
    """
    unique_seqs = list(dict.fromkeys(seq_list))
    if not unique_seqs:
        return {}
    return dict(_homology_worker(seq) for seq in unique_seqs)

# Peptide functions adapted from Mordret et al., Mol.Cell, 2019 & Tsour et al., Nature, 2026 
# to identify amino acid substitutions in MQ DP search results
def prot_nterm(sequence):
    """
    Does the peptide originate at the protein's N-term?
    """
    for start in SA_search(sequence, W_aa_ambiguous, sa_ambiguous):
        s = int(start)
        if W_aa_ambiguous[s-1] == '*' or W_aa_ambiguous[s-2] == '*':
            return True
    return False

def prot_cterm(sequence):
    """
    Does the peptide end at the protein's C-term?
    """
    l = len(sequence)
    for start in SA_search(sequence, W_aa_ambiguous, sa_ambiguous):
        if W_aa_ambiguous[int(start)+l] == '*':
            return True
    return False

def pep_cterm(ms):
    """
    Returns the probability that C term AA was modified.
    """
    if ms[-1] == ')':
        return float(ms[:-1].split('(')[-1]) >= 0.95
    return False

def pep_nterm(ms):
    """
    Returns the probability that N term AA was modified.
    """
    if len(ms) > 1 and ms[1] == '(':
        return float(ms[2:].split(')')[0]) >= 0.95
    return False

def refine_localization_probabilities(modified_seq):
    """
    Input: modified sequence (a string of AA with p of each to contain modification: APKV(.7)ML(.3)L means that V was modified with p = .7 and L with p = .3)
    Output: 2 lists: all candidate residues and their positional probabilities.
    """
    sites   = [modified_seq[m.start()-1] for m in re.finditer(r'\(', modified_seq)]
    weights = [float(x) for x in re.findall(r'\(([^\)]+)\)', modified_seq)]
    return sites, weights

def get_aa_subs(dp_dict, idx, tol=5, pp=0):
    """
    Find DP sequences with mass shift that is within tolerance of theoretical mass shift of AAS. Positional probability value must be reported by MQ for base residue consideration
    Input: idx = iterator value, tol = mass shift error tolerance, pp = positional probability threshold
    Output: AAS information
    """
    cand_res  = dp_dict['DP candidate residues'][idx]
    pos_probs = dp_dict['DP positional probabilities'][idx]
    DP_deltam = dp_dict['DP Mass Difference'][idx]
    DP_mz     = dp_dict['m/z'][idx]
    mtol      = DP_mz * (tol / 1e6)

    out_origin, out_idx, out_subs, out_prob, out_err = [], [], [], [], []
    for i, res in enumerate(cand_res):
        if pos_probs[i] <= pp:
            continue
        sub_keys  = subs_dict[res]
        masses    = subs_mass_dict[res]
        if len(masses) == 0:
            continue
        same_sign = ((DP_deltam > 0) == (masses > 0))
        in_tol    = (masses > DP_deltam - mtol) & (masses < DP_deltam + mtol)
        hits      = np.where(same_sign & in_tol)[0]
        for j in hits:
            out_origin.append(res)
            out_idx.append(i)
            out_subs.append(sub_keys[j])
            out_prob.append(pos_probs[i])
            out_err.append(float(abs(DP_deltam - masses[j])))
    return [out_origin, out_idx, out_subs, out_prob, out_err]

def get_mistranslated_seq(dp_dict, idx):
    """
    Input: idx = iterator value as loop over data dictionary
    Output: list of potential sequences that could explain mass shift
    """
    DP             = dp_dict['DP Probabilities'][idx]
    cand_res       = dp_dict['DP candidate residues'][idx]
    cand_res_idx   = [m.start()-1 for m in re.finditer(r'\(', DP)]
    sub_origin_idx = dp_dict['origin aa index'][idx]
    sub_dest       = dp_dict['destination aa'][idx]
    mps = []
    for i, res in enumerate(cand_res):
        if i in sub_origin_idx:
            j       = sub_origin_idx.index(i)
            ci      = cand_res_idx[i]
            mod_seq = re.sub(r'\([^\)]+\)', '', DP[:ci] + sub_dest[j] + DP[ci+1:])
            mps.append(mod_seq)
    return mps

def find_PTMs(dp_dict, idx, tol=10):
    """
    finds and annotates DPs that can be explained by PTM
    """
    cand_res  = dp_dict['DP candidate residues'][idx]
    pos_probs = dp_dict['DP positional probabilities'][idx]
    DP_deltam = dp_dict['DP Mass Difference'][idx]
    DP_mz     = dp_dict['m/z'][idx]
    pn        = dp_dict['Peptide N-term'][idx]
    pc        = dp_dict['Peptide C-term'][idx]
    mtol      = DP_mz * (tol / 1e6)
    PTM, PTM_aa, PTM_prob, PTM_dm = [], [], [], []
    for i, res in enumerate(cand_res):
        sub = mod_df.loc[(mod_df.site==res)|(mod_df.site=='N-term')|(mod_df.site=='C-term')]
        for _, row in sub.iterrows():
            dm = row['mass shift']
            if DP_deltam - mtol < dm < DP_deltam + mtol:
                pos = row['position']
                if (pos=='Anywhere') or (pos=='Any N-term' and pn) or (pos=='Any C-term' and pc):
                    PTM.append(row['modification']); PTM_aa.append(res)
                    PTM_prob.append(pos_probs[i]);    PTM_dm.append(DP_deltam-dm)
    return [PTM, PTM_aa, PTM_prob, PTM_dm]

### Main Processing 
print(f"Reading allPeptides.txt...")
_ALL_PEP_DF = pd.read_csv(
    f"{MQ_dir}allPeptides.txt", sep='\t', low_memory=False
)
print(f"allPeptides.txt loaded: {len(_ALL_PEP_DF)} rows.")

def get_data_dict(s):
    """ Function to identify and annotate peptides with AAS and other canonical PTMs """
    print(f"Processing sample: {s}")

    sample_df = _ALL_PEP_DF  
    if 'DP Base Sequence' not in sample_df.columns:
        raise KeyError(f"Missing 'DP Base Sequence'. Columns: {list(sample_df.columns)}")

    allPep_count = len(sample_df)
    dp_df = sample_df.loc[
        sample_df['DP Mass Difference'].notna() & (sample_df['DP PEP'] <= 0.01),
        ['Raw file','Charge','m/z','Mass','Mass precision [ppm]','Retention time',
         'Sequence','Proteins','Intensity','DP Mass Difference','DP Time Difference',
         'DP PEP','DP Base Sequence','DP Probabilities','DP Positional Probability',
         'DP Base Raw File','DP Base Scan Number','DP Mod Scan Number',
         'DP Proteins','DP Ratio mod/base']
    ].reset_index(drop=True)

    dp_dict = dp_df.to_dict()
    dp_dict['count allPeptides'] = allPep_count
    dp_dict['count DP']          = len(dp_df)

    # Identify candidate residues & modification probabilities
    crd = {i: refine_localization_probabilities(v)
           for i, v in dp_dict['DP Probabilities'].items()}
    dp_dict['DP candidate residues']         = {i: crd[i][0] for i in crd}
    dp_dict['DP positional probabilities']   = {i: crd[i][1] for i in crd}
    dp_dict['count candidate residues per peptide'] = {i: len(v) for i, v in dp_dict['DP candidate residues'].items()}

    dp_dict['Protein N-term'] = {i: prot_nterm(v) for i, v in dp_dict['DP Base Sequence'].items()}
    dp_dict['Protein C-term'] = {i: prot_cterm(v) for i, v in dp_dict['DP Base Sequence'].items()}
    dp_dict['Peptide N-term'] = {i: pep_nterm(v)  for i, v in dp_dict['DP Probabilities'].items()}
    dp_dict['Peptide C-term'] = {i: pep_cterm(v)  for i, v in dp_dict['DP Probabilities'].items()}

    # Annotate DPs that can be explained by AAS 
    asd = {i: get_aa_subs(dp_dict, i) for i in dp_dict['DP candidate residues']}
    dp_dict['origin aa']                    = {i: asd[i][0] for i in asd}
    dp_dict['origin aa index']              = {i: asd[i][1] for i in asd}
    dp_dict['aa subs']                      = {i: asd[i][2] for i in asd}
    dp_dict['aa subs positional probability']= {i: asd[i][3] for i in asd}
    dp_dict['aa subs mass error (ppm)']     = {i: asd[i][4] for i in asd}
    dp_dict['destination aa']               = {i: [x[-1] for x in v] if v else []
                                               for i, v in dp_dict['aa subs'].items()}
    dp_dict['mistranslated sequence']       = {i: get_mistranslated_seq(dp_dict, i)
                                               for i in dp_dict['aa subs']}
    dp_dict['count mistranslated peptides'] = sum(1 for x in dp_dict['origin aa'].values() if x)
    dp_dict['count aa subs per peptide']    = {i: len(v) for i, v in dp_dict['origin aa'].items()}

    # Annotate DPs that can be explained by canonical PTM
    ptmd = {i: find_PTMs(dp_dict, i) for i in dp_dict['DP candidate residues']}
    dp_dict['PTM']                                      = {i: ptmd[i][0] for i in ptmd}
    dp_dict['PTM site']                                 = {i: ptmd[i][1] for i in ptmd}
    dp_dict['PTM positional probability']               = {i: ptmd[i][2] for i in ptmd}
    dp_dict['PTM mass error [observed-expected] (ppm)'] = {i: ptmd[i][3] for i in ptmd}

    # Annotate whether or not peptide with AAS has homologous sequence that could have arisen from elsewhere in genome
    all_seqs = [x for v in dp_dict['mistranslated sequence'].values() for x in v]
    print(f"{s}: {len(all_seqs)} MTP sequences to check for homology "
                 f"({len(set(all_seqs))} unique).")

    homology_cache = batch_homology_search(all_seqs)   # {seq: bool}

    # Map cached results back to dp_dict entries
    dp_dict['genome_homolog'] = {
        i: [homology_cache.get(x, False) for x in v]
        for i, v in dp_dict['mistranslated sequence'].items()
    }
    
    for frame_col in [f'{f}-frame genome substring' for f in range(1,7)]:
        dp_dict[frame_col] = dp_dict['genome_homolog']

    print(f"Processing for sample {s} completed.")
    return dp_dict

### Run pipeline
if __name__ == '__main__':
    fout = open('aas_detect.out', 'w')
    fout.write('starting AAS pipeline\n')

    # map over samples
    with mp.Pool() as pool:
        data_dict_list = pool.map(get_data_dict, tissues)
    data_dict = {s: data_dict_list[i] for i, s in enumerate(tissues)}
    pickle.dump(data_dict, open(aas_dir+'DP_dict.p', 'wb'))

    fout.write('filtering data dict\n')
    mtp_dict = {}
    ptm_dict = {}

    for s in tissues:
        dp_dict = data_dict[s]
        fout.write(s + '\n')
        ptm_idx = [i for i in dp_dict['DP PEP'] if len(dp_dict['PTM'][i]) > 0]
        mtp_idx = [i for i in dp_dict['DP PEP']
                   if not dp_dict['PTM'][i] and dp_dict['aa subs'][i]]
        curr_mtp, curr_ptm = {}, {}
        for k, v in dp_dict.items():
            if k not in ('count mistranslated peptides','count allPeptides','count DP'):
                curr_mtp[k] = {i: v[i] for i in mtp_idx}
                curr_ptm[k] = {i: v[i] for i in ptm_idx}
        mtp_dict[s] = curr_mtp
        ptm_dict[s] = curr_ptm
        print(f'{s}: N MTPs = {len(mtp_dict[s]["mistranslated sequence"])}')

    pickle.dump(mtp_dict, open(aas_dir+'MTP_dict.p', 'wb'))
    pickle.dump(ptm_dict, open(aas_dir+'PTM_dict.p', 'wb'))
    
    fout.write('files saved successfully!\n')
    fout.close()
    print("AAS pipeline completed!")