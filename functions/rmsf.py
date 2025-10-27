# Importing the necessary libraries for plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 

from moleculekit.molecule import Molecule 
import moleculekit.projections.metricfluctuation as fluct
import sys
import os
import re

def rmsf(filename=str, topology=str, trajectory=str, atomsel=str, refmol=None, trajalnsel=str, refalnsel='protein and name CA'):
    """
    Compute and plot the Root Mean Square Fluctuation (RMSF) for selected atoms in a trajectory.

    This function uses MoleculeKit's `MetricFluctuation` to compute RMSF values for 
    a given atom selection, relative to either a reference molecule or the trajectory mean.
    Results are stored in a CSV file for reuse.

    Parameters
    ----------
        filename: str
            Prefix for saving results (CSV and PDF).
        topology: str 
            Path to the topology file (e.g., PDB, PSF).
        trajectory: str 
            Path to the trajectory file (e.g., DCD, XTC).
        atomsel: str 
            Atom selection string for which RMSF will be computed 
            (e.g., `"protein and name CA"` for alpha carbons).
        refmol: str, optional 
            Path to a reference structure file. 
            If `None`, the trajectory mean will be used as the reference.
        trajalnsel: str Atom selection string for aligning the trajectory.
        refalnsel: str, optional 
            Atom selection string for aligning the reference. 
            Defaults to `'protein and name CA'`.

    Output
    -------
    - CSV file: `../results/tables/{filename}_rmsf.csv`
    """
    csv_file = f"../results/tables/{filename}_rmsf.csv"
    
    if os.path.exists(csv_file):
        # If CSV exists, just load it
        data = pd.read_csv(csv_file)
    else:
        if refmol == None:
            ref_mol = None
            print('RMSF will be computed around the trajectory mean atomic positions, else add a refmol Molecule object.')
        else:
            ref_mol = Molecule(refmol, validateElements = False)
            #ref_mol.filter("not resname TIP3") #run filter_prod.py before calling this function
        mol = Molecule(topology, validateElements = False)
        mol.read(trajectory)
        mol.filter("not resname TIP3")
        #mol.dropFrames(keep = np.arange(5000, 10000)) #keep only the second half
    

        #compute rmsf
        #superpose muovere tutta proteina su refalnsel
        mol.align(sel = refalnsel)
        met=fluct.MetricFluctuation(atomsel = atomsel, refmol = ref_mol, trajalnsel= trajalnsel, refalnsel = refalnsel, mode = 'atom', pbc =False) #allowed also "residue" 
        fluct_values = met.project(mol)

        rmsf = np.sqrt(fluct_values.mean(axis=0))

        #residue index (assume 1 atom per residue i.e. CA)
        resid= mol.get('resid', atomsel) #its a dimer, numbers will be repeated.
        segid= mol.get('segid', atomsel)
        #store data in csv
        data = pd.DataFrame({
            'resid': resid,
            'segid': segid,
            'rmsf': rmsf
        })
        data.to_csv(csv_file, index=False)

    return
