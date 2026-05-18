#!/bin/usr python

import os
import pickle
import warnings
import logging
import time
import numpy as np
from itertools import groupby
from operator import itemgetter
from Bio import SeqIO
from Bio import BiopythonWarning

warnings.simplefilter('ignore',BiopythonWarning)

# Build genome frame translations for decode pipeline
# Script based on frame translation generation in Tsour et al. Nature, 2026
# Last updated by: Alex Maropakis, 05-18-2026

print("===================================================")
print("Building genome frame translations")
print("===================================================")

scratch_dir='/scratch/maropakis.a/'

GENOMES={
'human':{
'fasta':scratch_dir+'Dependencies/FASTA/HUMAN_GENOME.fna',
'outdir':scratch_dir+'Dependencies/frame_translations/human/'
},
'mouse':{
'fasta':scratch_dir+'Dependencies/FASTA/MOUSE_GENOME.fna',
'outdir':scratch_dir+'Dependencies/frame_translations/mouse/'
}
}

###############################################################################
# Logging
###############################################################################

log_file=scratch_dir+'Dependencies/frame_translations/gen_frameshifts.log'

logging.basicConfig(filename=log_file,level=logging.INFO,format='%(asctime)s | %(levelname)s | %(message)s')

logging.info('Starting frame translation generation')

###############################################################################
# Suffix array
###############################################################################

def suffix_array(text,_step=16):
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
    >>> suffix_array(text='banana')
    ([5, 3, 1, 0, 4, 2], [3, 2, 5, 1, 4, 0], [0, 1, 3, 0, 0, 2])
    
    Explanation: 'a' < 'ana' < 'anana' < 'banana' < 'na' < 'nana'
    The Longest Common String is 'ana': lcp[2] == 3 == len('ana')
    It is between  tx[sa[1]:] == 'ana' < 'anana' == tx[sa[2]:]
    """

    tx=text
    size=len(tx)
    step=min(max(_step,1),len(tx))
    sa=list(range(len(tx)))
    sa.sort(key=lambda i:tx[i:i+step])
    grpstart=size*[False]+[True]
    rsa=size*[None]
    stgrp,igrp='',0

    for i,pos in enumerate(sa):

        st=tx[pos:pos+step]

        if st!=stgrp:
            grpstart[igrp]=(igrp<i-1)
            stgrp=st
            igrp=i

        rsa[pos]=igrp
        sa[i]=pos

    grpstart[igrp]=(igrp<size-1 or size==0)

    while grpstart.index(True)<size:

        nextgr=grpstart.index(True)

        while nextgr<size:

            igrp=nextgr
            nextgr=grpstart.index(True,igrp+1)
            glist=[]

            for ig in range(igrp,nextgr):

                pos=sa[ig]

                if rsa[pos]!=igrp:
                    break

                newgr=rsa[pos+step] if pos+step<size else -1
                glist.append((newgr,pos))

            glist.sort()

            for ig,g in groupby(glist,key=itemgetter(0)):

                g=[x[1] for x in g]
                sa[igrp:igrp+len(g)]=g
                grpstart[igrp]=(len(g)>1)

                for pos in g:
                    rsa[pos]=igrp

                igrp+=len(g)

        step*=2

    del grpstart
    del rsa

    return sa

###############################################################################
# Build translations
###############################################################################

def build_ref_translation(path_to_fasta,frameshift_dir,f=[1,2,3,4,5,6]):
    """ build reference fasta for removal of homologous sequences """
    """ f = frame """

    logging.info(f'Starting frame {f}')

    record_list=[]
    translated_record_list=[]
    record_dict={}
    boundaries_aa=[0]

    t0=time.time()

    for i,record in enumerate(SeqIO.parse(open(path_to_fasta,'r+'),'fasta')):

        if i%100==0:
            print(f'Processing record {i}')
            logging.info(f'Frame {f} | record {i}')

        if f in [1,2,3]:
            record.seq=record.seq.upper()[f:]
        else:
            record.seq=record.seq.upper()[::-1][f:]

        translation=str(record.seq.translate())
        bits=record.description.split(' ')

        record_list.append(record)
        translated_record_list.append(translation)

        boundaries_aa.append(boundaries_aa[-1]+len(translation))

    boundaries_aa=np.array(boundaries_aa[1:])

    print("Concatenating translations...")
    logging.info(f'Frame {f} concatenating translations')

    W_aa=''.join(translated_record_list)

    print(f'Translated AA length: {len(W_aa):,}')
    logging.info(f'Frame {f} translated AA length={len(W_aa)}')

    print("Building suffix array...")
    logging.info(f'Frame {f} building suffix array')

    sa=suffix_array(W_aa)

    print("Replacing I -> L...")
    logging.info(f'Frame {f} replacing I/L ambiguity')

    W_aa_ambiguous=W_aa.replace('I','L')
    sa_ambiguous=suffix_array(W_aa_ambiguous)

    print("Saving pickles...")
    logging.info(f'Frame {f} saving pickles')

    pickle.dump(W_aa_ambiguous,open(frameshift_dir+'W'+str(f)+'_aa_ambig.p','wb'))
    pickle.dump(sa_ambiguous,open(frameshift_dir+'s'+str(f)+'a_ambig.p','wb'))

    runtime=(time.time()-t0)/60

    print(f'Finished frame {f}')
    logging.info(f'Finished frame {f} | runtime={runtime:.2f} min')

    return(W_aa_ambiguous,sa_ambiguous)

###############################################################################
# Main
###############################################################################

for species,cfg in GENOMES.items():

    fasta_path=cfg['fasta']
    outdir=cfg['outdir']

    print("===================================================")
    print(f"Species: {species}")
    print(f"FASTA:   {fasta_path}")
    print(f"Output:  {outdir}")
    print("===================================================")

    logging.info(f'Species={species}')
    logging.info(f'FASTA={fasta_path}')
    logging.info(f'Output={outdir}')

    os.makedirs(outdir,exist_ok=True)

    for frame in range(1,7):

        wp=os.path.join(outdir,f'W{frame}_aa_ambig.p')
        sp=os.path.join(outdir,f's{frame}a_ambig.p')

        if os.path.exists(wp) and os.path.exists(sp):

            print(f'Frame {frame} already exists. Skipping.')
            logging.info(f'Skipping existing frame {frame}')

            continue

        build_ref_translation(fasta_path,outdir,frame)

print("===================================================")
print("All frame translations completed.")
print("===================================================")

logging.info('All frame translations completed')