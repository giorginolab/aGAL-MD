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
    
    csv_file = f"../results/{filename}_rmsf.csv"
    
    if os.path.exists(csv_file):
        # If CSV exists, just load it
        data = pd.read_csv(csv_file)
    else:
        if refmol == None:
            print('RMSF will be computed around the trajectory mean, Else add a refmol Molecule object.')
        ref_mol = Molecule(refmol, validateElements = False)
        ref_mol.filter("not resname TIP3")
        mol = Molecule(topology, validateElements = False)
        mol.read(trajectory)
        mol.filter("not resname TIP3")
        time_values = np.arange(0,len(mol.time))

        #compute rmsf
        met=fluct.MetricFluctuation(atomsel = atomsel, refmol = ref_mol, trajalnsel= trajalnsel, refalnsel = refalnsel, mode = 'atom') #allowed also "residue" 
        fluct_values = met.project(mol)

        rmsf = np.sqrt(fluct_values.mean(axis=0))

        #residue index (assume 1 atom per residue i.e. CA)
        resid= mol.get('resid', atomsel) #its a dimer, numbers will be repeated.
        chain = mol.get('chain', atomsel)
        #store data in csv
        data = pd.DataFrame({
            'resid': resid,
            'chain': chain,
            'rmsf': rmsf
        })
        data.to_csv(csv_file, index=False)

    #plotting
    fig,ax = plt.subplots(figsize=(10,5))

    ax.plot(data['resid'],data['rmsf'], color="black", linestyle="-")

    ax.set_xlabel("Sequence")
    ax.set_ylabel(r"C$_\alpha$ RMSF (nm)")

    fig.savefig(f"../results/{filename}_rmsf.png", dpi=300)
    plt.show()

    #just in case
    return data