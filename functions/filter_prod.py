#!/usr/bin/env python3

#REMOVE WATER AND SLICE TRAJ EVERY 10th FRAME

#outside the protein folder (i.e. within 3GXT_reglyco/)
#python ../functions/filter_prod.py
import sys
import os
from moleculekit.molecule import Molecule

replica_folders = sorted([d for d in os.listdir('.') if os.path.isdir(d)])

if not replica_folders: 
    print("No replica folders found in current directory.")
    sys.exit(1)

for replica in replica_folders:
    target_folder = os.path.join(replica, 'production')
    
    psf_path = os.path.join(target_folder, 'structure.psf')
    xtc_path = os.path.join(target_folder, 'output.xtc')
    output_psf = os.path.join(replica, f'{replica}_noh.psf')
    output_xtc = os.path.join(replica, f'{replica}_noh_s10.xtc')
        
    if not os.path.isdir(target_folder):
        print(f'Missing production folder in {replica} directory.')

    if os.path.exists(output_psf) and os.path.exists(output_xtc):
        print(f"Skipping {replica}: filtered files already exist.")
        continue

    #works after mol_prep.ipynb
    mol = Molecule(psf_path, validateElements=False)
    mol.read(xtc_path, skip=10) #shorten traj
    mol.filter('not resname TIP3')  #remove water

    mol.write(output_psf)
    mol.write(output_xtc)

    print(f'Completed filtering for molecule {replica}.')