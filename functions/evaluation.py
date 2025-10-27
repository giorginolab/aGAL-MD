#!/usr/bin/env python

import pandas as pd 
import numpy as np
from moleculekit.molecule import Molecule

from rmsd import * 
from rmsf import *


#############

'''
# RMSD - WHOLE PROTEIN SPLIT BY MONOMER
structure = ['apo', 'apo_N215S', 'apo_R301Q', 'DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh_s10.xtc'
        rmsd(f'{s}{r}_CA_P0', topology, trajectory,  trajrmsdstr='protein and name CA and segid P0', trajalnstr='protein and name CA',refalnstr='protein and name CA',  refrmsdstr = 'protein and name CA and segid P0')
        rmsd(f'{s}{r}_CA_P1', topology, trajectory,  trajrmsdstr='protein and name CA and segid P1', trajalnstr='protein and name CA',refalnstr='protein and name CA', refrmsdstr = 'protein and name CA and segid P1') 

'''
     

#############


'''
#RMSD SEPARATE LIGAND
structure = ['DGJ', 'DGJ_N215S', 'DGJ_R301Q'] #apo doesn't have it obviously
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh.psf' #non caricare come topology COME PSF
        trajectory = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh_s10.xtc'
        
        rmsd(f'{s}{r}_lig_P0', topology, trajectory,  trajrmsdstr='resname DGJ and resid 1 and noh',trajalnstr='segid P0 and name CA',refalnstr='segid P0 and name CA', refrmsdstr = 'resname DGJ and resid 1 and noh' )
        rmsd(f'{s}{r}_lig_P1', topology, trajectory,  trajrmsdstr='resname DGJ and resid 2 and noh',trajalnstr='segid P1 and name CA',refalnstr='segid P1 and name CA', refrmsdstr = 'resname DGJ and resid 2 and noh' )
'''

#############

################################## RMSF ###################################

#############


# RMSF - protein monomer
'''structure = ['apo', 'apo_N215S', 'apo_R301Q', 'DGJ', 'DGJ_N215S', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh_s10.xtc'
        
        rmsf(f'{s}{r}_CA_P0', topology, trajectory, atomsel = 'segid P0 and name CA', trajalnsel='segid P0 and name CA', refalnsel='segid P0 and name CA' )
        rmsf(f'{s}{r}_CA_P1', topology, trajectory, atomsel = 'segid P1 and name CA', trajalnsel='segid P1 and name CA', refalnsel='segid P1 and name CA' )
'''

#############
'''
# RMSF - GLYCANS IN RESPECT TO PROTEIN
structure = ['apo', 'apo_R301Q', 'DGJ', 'DGJ_R301Q']
replica = ['_1', '_2', '_3']

for s in structure:
    for r in replica:
        topology = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh.psf'
        trajectory = f'../3GXT_reglyco/{s}{r}/{s}{r}_noh_s10.xtc'
        
        rmsf(f'{s}{r}_gly_P2', topology, trajectory, atomsel = 'segid P2', trajalnsel='segid P0 and name CA', refalnsel='segid P0 and name CA' )
        rmsf(f'{s}{r}_gly_P3', topology, trajectory, atomsel = 'segid P3', trajalnsel='segid P0 and name CA', refalnsel='segid P0 and name CA' )
        rmsf(f'{s}{r}_gly_P4', topology, trajectory, atomsel = 'segid P4', trajalnsel='segid P0 and name CA', refalnsel='segid P0 and name CA' )
        rmsf(f'{s}{r}_gly_P5', topology, trajectory, atomsel = 'segid P5', trajalnsel='segid P1 and name CA', refalnsel='segid P1 and name CA' )
        rmsf(f'{s}{r}_gly_P6', topology, trajectory, atomsel = 'segid P6', trajalnsel='segid P1 and name CA', refalnsel='segid P1 and name CA' )
        rmsf(f'{s}{r}_gly_P7', topology, trajectory, atomsel = 'segid P7', trajalnsel='segid P1 and name CA', refalnsel='segid P1 and name CA' )


mutant=['apo_N215S', 'DGJ_N215S']
for m in mutant:
    for r in replica:
        topology = f'../3GXT_reglyco/{m}{r}/{s}{r}_noh.psf'
        trajectory = f'../3GXT_reglyco/{m}{r}/{s}{r}_noh_s10.xtc'
        rmsf(f'{m}{r}_gly_P2', topology, trajectory, atomsel = 'segid P2', trajalnsel='segid P0 and name CA', refalnsel='segid P0 and name CA' )
        rmsf(f'{m}{r}_gly_P3', topology, trajectory, atomsel = 'segid P3', trajalnsel='segid P0 and name CA', refalnsel='segid P0 and name CA' )
        
        rmsf(f'{m}{r}_gly_P5', topology, trajectory, atomsel = 'segid P5', trajalnsel='segid P1 and name CA', refalnsel='segid P1 and name CA' )
        rmsf(f'{m}{r}_gly_P6', topology, trajectory, atomsel = 'segid P6', trajalnsel='segid P1 and name CA', refalnsel='segid P1 and name CA' )
'''