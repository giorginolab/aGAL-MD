#!/usr/bin/env python

import pandas as pd 
import numpy as np
from moleculekit.molecule import Molecule

from rmsd import * 
from rmsf import *


'''
# RMSD - WHOLE PROTEIN - CA
structure = ['apo', 'apo_N215S', 'apo_R301Q', 'DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/production/structure.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/production/output.xtc'
        
        rmsd(f'{s}{r}_CA', topology, trajectory, trajrmsdstr='protein and name CA',trajalnstr='protein and name CA',refalnstr='protein and name CA', refrmsdstr = 'protein and name CA', title = True)
'''       
'''
# RMSD - WHOLE PROTEIN SPLIT BY MONOMER
structure = ['apo', 'apo_N215S', 'apo_R301Q', 'DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/production/structure.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/production/output.xtc'
        rmsd(f'{s}{r}_CA_chain_A', topology, trajectory,  trajrmsdstr='protein and name CA and chain A', refrmsdstr = 'protein and name CA and chain A')
        rmsd(f'{s}{r}_CA_chain_B', topology, trajectory,  trajrmsdstr='protein and name CA and chain B', refrmsdstr = 'protein and name CA and chain B')
 
'''
 
# RMSF - WHOLE PROTEIN  
structure = ['apo', 'apo_N215S', 'apo_R301Q', 'DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/production/structure.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/production/output.xtc'
        
        rmsf(f'{s}{r}_CA', topology, trajectory, atomsel = 'protein and name CA', trajalnsel='protein and name CA', refalnsel='protein and name CA', title = True)
         

'''
# RMSD - LIGAND IN RESPECT TO LIGAND 

structure = ['DGJ', 'DGJ_N215S', 'DGJ_R301Q'] #apo doesn't have it obviously
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/production/structure.psf' #non caricare come topology COME PSF
        trajectory = f'../3GXT_reglyco/{s}{r}/production/output.xtc'
        
        rmsd(f'{s}{r}_lig', topology, trajectory,  trajrmsdstr='resname DGJ',trajalnstr='protein and name CA',refalnstr='protein and name CA', refrmsdstr = 'resname DGJ', title = True)
'''


# RMSF -LIGAND IN RESPECT TO PROTEIN

structure = ['DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/production/structure.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/production/output.xtc'
        
        rmsf(f'{s}{r}_lig', topology, trajectory, atomsel = 'resname DGJ', trajalnsel='protein and name CA', refalnsel='protein and name CA', title = True)   


# RMSF - GLYCANS IN RESPECT TO PROTEIN
structure = ['apo', 'apo_N215S', 'apo_R301Q', 'DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/production/structure.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/production/output.xtc'
        
        rmsf(f'{s}{r}_gly', topology, trajectory, atomsel = 'resname AMAN or resname BGLC', trajalnsel='protein and name CA', refalnsel='protein and name CA', title = True)
