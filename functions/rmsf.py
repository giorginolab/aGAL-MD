# Importing the necessary libraries for plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 

from moleculekit.molecule import Molecule 
import moleculekit.projections.metricfluctuation as fluct
import sys
import os
import re

def rmsf(filename=str, topology=str, trajectory=str, atomsel=str, refmol=None, trajalnsel=str, refalnsel='protein and name CA', title = False):
    """
    Compute and plot the Root Mean Square Fluctuation (RMSF) for selected atoms in a trajectory.

    This function uses MoleculeKit's `MetricFluctuation` to compute RMSF values for 
    a given atom selection, relative to either a reference molecule or the trajectory mean.
    Results are stored in a CSV file for reuse, and a PDF plot of RMSF values per residue 
    is generated and saved.

    Args:
        filename (str): Prefix for saving results (CSV and PDF).
        topology (str): Path to the topology file (e.g., PDB, PSF).
        trajectory (str): Path to the trajectory file (e.g., DCD, XTC).
        atomsel (str): Atom selection string for which RMSF will be computed 
            (e.g., `"protein and name CA"` for alpha carbons).
        refmol (str or None, optional): Path to a reference structure file. 
            If `None`, the trajectory mean will be used as the reference.
        trajalnsel (str): Atom selection string for aligning the trajectory.
        refalnsel (str, optional): Atom selection string for aligning the reference. 
            Defaults to `'protein and name CA'`.
        title (bool, optional): If True, the plot will include the filename as the title. 
            Defaults to False.

    Returns:
        None
            Saves the following outputs:
            - `../results/{filename}_rmsf.csv`: Table with `resid`, `chain`, and `rmsf`.
            - `../results/{filename}_rmsf.pdf`: RMSF plot comparing chains A and B.

    Notes:
        - Glycosylation sites (residues 139, 192, 215) are annotated as vertical lines.
        - Assumes TIP3 (water) residues are filtered out.
        - Assumes one atom per residue in `atomsel` (e.g., CA atoms) for indexing.
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
            ref_mol.filter("not resname TIP3")
        mol = Molecule(topology, validateElements = False)
        mol.read(trajectory)
        mol.filter("not resname TIP3")
        mol.dropFrames(keep = np.arange(5000, 10000)) #keep only the second half
    

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

    #data["resid_chain"] = data["resid"].astype(str) + data["chain"]

#plotting
    fig, ax = plt.subplots(figsize=(20,5))

    chainA = data[data['chain']=='A']
    chainB = data[data['chain']=='B']

    ax.plot(chainA['resid'], chainA['rmsf'], linestyle="-", label='Monomer A')
    ax.plot(chainB['resid'], chainB['rmsf'], linestyle="-", label='Monomer B')

    # Glycosylation sites
    first = True
    for marker_resid in [139, 192, 215]:
        ax.axvline(marker_resid, color="black", linestyle="-.", alpha=0.7, label='Glycosylation site' if first else "")
        ymax = ax.get_ylim()[1]
        ax.text(marker_resid + 0.5, 0.5, str(marker_resid), rotation=0,va='center', ha='left', fontsize=8, color='black')
        first = False

    # X-axis ticks
    step = 10
    unique_resid = np.unique(data['resid'])
    #special_resid = [139, 192, 215] #if you want to include add in the next line
    all_ticks = sorted(set(list(unique_resid[::step])))
    ax.set_xticks(all_ticks)
    ax.set_xticklabels(all_ticks, rotation=45, fontsize=8)
    ax.set_xlim(data['resid'].min() - 0.5, data['resid'].max() + 0.5)
    #ax.set_ylim(0,8)
    ax.set_xlabel("Residue ID")
    ax.set_ylabel("RMSF (Å)")
    plt.legend()

    if title:
        plt.title(filename)

    fig.savefig(f"../results/single_plots/{filename}_rmsf.pdf", dpi=300)
    plt.close()

    return