#!/usr/bin/env python3

"""
msconvert.py 

Convert every .raw file in raw_dir with ThermoRawFileParser (-f=2 mzML) into mzML, 
skipping any that already exist. Then symlink mzML files into spectra root. 

python msconvert.py /scratch/maropakis.a/MQ_raw/____ \
    --plex ____ \
    --spectra-root /scratch/maropakis.a/spectra

"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import Counter
import shutil

def convert_raws(raw_dir, trfp):
    """Convert every .raw in raw_dir to plain mzML in place; skip existing. Return mzML paths."""
    mzmls = []
    raws = sorted(glob.glob(os.path.join(raw_dir, '*.raw')))
    if not raws:
        # why: raw_dir may already hold mzML (pre-converted); fall through to collect those.
        mzmls = sorted(glob.glob(os.path.join(raw_dir, '*.mzML')))
        print(f'  no .raw in {raw_dir}; found {len(mzmls)} existing .mzML')
        return mzmls
    for raw in raws:
        base = os.path.splitext(os.path.basename(raw))[0]
        out = os.path.join(raw_dir, f'{base}.mzML')
        if os.path.getsize(out) if os.path.exists(out) else 0:
            print(f'  SKIP convert {base} (mzML exists)')
        else:
            # why: -f=2 = plain indexed mzML; -g would gzip and FragPipe silently skips .mzML.gz.
            subprocess.run([trfp, f'-i={raw}', f'-o={raw_dir}', '-f=2', '-l=3'], check=True)
            print(f'  converted {base}')
        if os.path.exists(out):
            mzmls.append(out)
    return mzmls

def stage_spectra(mzmls, spectra_root, plex):
    """Symlink mzML into <spectra-root>/<plex>/. Return staged dir."""
    dst = os.path.join(spectra_root, plex)
    os.makedirs(dst, exist_ok=True)
    for f in mzmls:
        link = os.path.join(dst, os.path.basename(f))
        if not os.path.lexists(link):
            os.symlink(os.path.abspath(f), link)
    print('staged mzML')
    return dst 

def main():
    ap = argparse.ArgumentParser(description='Convert .raw into .mzML')
    ap.add_argument('raw_dir', help='dir holding plex .raw files')
    ap.add_argument('--plex', required=True, help='plex token, e.g. acgb1 / pooled / aorta')
    ap.add_argument('--spectra-root', required=True)
    ap.add_argument('--trfp', default=os.path.expanduser('~/thermoRawFileParser/ThermoRawFileParser'))
    a = ap.parse_args()

    if a.raw_dir:
        mzmls = convert_raws(a.raw_dir, a.trfp)
    if not mzmls:
        sys.exit('no mzML to stage.')
    
    stage_spectra(mzmls, a.spectra_root, a.plex)
    staged = sorted(glob.glob(os.path.join(a.spectra_root, a.plex, '*.mzML')))

if __name__ == '__main__':
    main()