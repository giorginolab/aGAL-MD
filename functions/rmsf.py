import os
import numpy as np
import pandas as pd
from moleculekit.molecule import Molecule
import moleculekit.projections.metricfluctuation as fluct


def rmsf(
    filename,
    topology,
    trajectory,
    atomsel,
    refmol=None,
    trajalnsel="protein and name CA",
    refalnsel="protein and name CA",
):
    """
    Compute the root-mean-square fluctuation (RMSF) of a trajectory with 
    respect to:
    - refmol, a provided reference structure;
    - the mean atomic position if no reference is provided.
    Results saved as CSV tables.
    
    If called twice with the same filename it will not recompute the RMSF
    if the corresponding CSV file is already present.

    This function uses MoleculeKit's `MetricFluctuation` to compute RMSF.

    Parameters
    ----------
    filename : str
        Base name for the output file stored under `../results/tables`.
    topology : str
        Path to the topology file (e.g., PDB, PSF).
    trajectory : str
        Path to the trajectory file (e.g., DCD, XTC).
    atomsel : str
        Atom selection string for which RMSF will be computed
        (e.g., "protein and name CA" for alpha carbons).
    refmol : str or None, optional
        Path to a reference structure file.
        If `None`, the trajectory mean will be used as the reference.
    trajalnsel : str, optional
        Atom selection string for aligning the trajectory.
        Defaults to "protein and name CA".
    refalnsel : str, optional
        Atom selection string for aligning the reference.
        Defaults to "protein and name CA".

    Output
    -------
    - CSV file: `../results/tables/{filename}_rmsf.csv`
    """
    csv_file = f"../results/tables/{filename}_rmsf.csv"
    
    if os.path.exists(csv_file):
        # If CSV exists, just load it to confirm readability and exit.
        pd.read_csv(csv_file)
    else:
        if refmol is None:
            ref_mol = None
            print('RMSF will be computed around the trajectory mean atomic positions, else add a refmol Molecule object.')
        else:
            print('If not removed yet, run "filter_prod.py" to remove water before loading.')
            ref_mol = Molecule(refmol, validateElements=False)

        mol = Molecule(topology, validateElements=False)
        mol.read(trajectory)
        mol.filter("not resname TIP3")

        # Compute RMSF
        mol.align(sel=refalnsel)
        met = fluct.MetricFluctuation(
            atomsel=atomsel,
            refmol=ref_mol,
            trajalnsel=trajalnsel,
            refalnsel=refalnsel,
            mode="atom",  # "residue" is also allowed if needed
            pbc=False,
        )
        fluct_values = met.project(mol)

        rmsf = np.sqrt(fluct_values.mean(axis=0))

        # Residue index (assume one atom per residue, i.e., CA)
        resid = mol.get("resid", atomsel)
        segid = mol.get("segid", atomsel)
        
        # Store data in CSV
        data = pd.DataFrame({
            'resid': resid,
            'segid': segid,
            'rmsf': rmsf
        })
        data.to_csv(csv_file, index=False)

    return
