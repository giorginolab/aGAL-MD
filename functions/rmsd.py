# Importing the necessary libraries for plotting
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd 

from moleculekit.molecule import Molecule 
import moleculekit.projections.metricrmsd as metricrmsd
import sys
import os
import re

def rmsd(filename=str, topology=str, trajectory=str, refmol=str, trajrmsdstr=str, refrmsdstr='protein and name CA'):
    """
    Compute and plot the root-mean-square deviation (RMSD) of a trajectory with 
    respect to a reference structure, saving results to CSV and PNG files.

    This function calculates the RMSD of a trajectory relative to a reference
    molecule using the MoleculeKit `MetricRmsd` projection. Results are saved in
    the `../results/` directory. If the CSV file already exists, the function 
    skips the computation and directly loads the saved data for plotting.

    Parameters
    ----------
    filename : str
        Base name for the output files (CSV and PNG) stored under `../results/`.
    topology : str
        Path to the topology file (e.g., PDB or PSF).
    trajectory : str
        Path to the trajectory file (e.g., DCD).
    refmol : str
        Path to the reference structure file used for RMSD calculation.
    trajrmsdstr : str
        Atom selection string for trajectory atoms used in RMSD calculation.
    refrmsdstr : str, optional
        Atom selection string for reference structure atoms (default: 'protein and name CA').

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns:
        - 'time' : trajectory time steps (ps)
        - 'rmsd' : computed RMSD values (nm)

    Outputs
    -------
    - CSV file: `../results/{filename}_rmsd.csv`
    - PNG file: `../results/{filename}_rmsd.png`

    Notes
    -----
    - Water molecules with resname `TIP3` are filtered out before calculations.
    - If the CSV file already exists, RMSD values are not recomputed.
    - The plot shows RMSD as a function of simulation time.
    """
    csv_file = f"../results/{filename}_rmsd.csv"
    
    if os.path.exists(csv_file):
        # If CSV exists, just load it
        data = pd.read_csv(csv_file)
    else:
        ref_mol = Molecule(refmol, validateElements = False)
        ref_mol.filter("not resname TIP3")
        mol = Molecule(topology, validateElements = False)
        mol.read(trajectory)
        mol.filter("not resname TIP3")
        time_values = np.arange(0,len(mol.time))

        #compute rmsd
        met=metricrmsd.MetricRmsd(ref_mol, trajrmsdstr= trajrmsdstr, refrmsdstr = refrmsdstr, refalnstr=refrmsdstr) 
        rmsd_values = met.project(mol)

        #store data in csv
        data = pd.DataFrame({
            'time': time_values,
            'rmsd': rmsd_values
        })
        data.to_csv(csv_file, index=False)

    #plotting
    fig, ax = plt.subplots(figsize=(10,5))

    ax.plot(data['time'],data['rmsd'], color="black", linestyle="-")

    ax.set_xlabel("time $t$ (ps)")
    ax.set_ylabel(r"C$_\alpha$ RMSD (nm)")

    fig.savefig(f"../results/{filename}_rmsd.png", dpi=300)
    plt.show()

    #just in case
    return data