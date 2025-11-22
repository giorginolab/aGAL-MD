import os
import numpy as np
import pandas as pd
from moleculekit.molecule import Molecule
import moleculekit.projections.metricrmsd as metricrmsd


def rmsd(
    filename,
    topology,
    trajectory,
    trajrmsdstr,
    trajalnstr=None,
    refrmsdstr="protein and name CA",
    refalnstr="protein and name CA",
    centerstr="protein and name CA",
):
    """
    Compute the root-mean-square deviation (RMSD) of a trajectory with 
    respect to a reference structure obtained from the frame 0 of the loaded trajectory, saving results to CSV tables.
    
    If called twice with the same filename it will not recompute the RMSD
    if the corresponding CSV file is already present.

    This function uses MoleculeKit's projection method `MetricRmsd` to compute RMSD.

    Parameters
    ----------
    filename : str
        Base name for the output file stored under `../results/tables`.
    topology : str
        Path to the topology file (e.g., PDB or PSF).
    trajectory : str
        Path to the trajectory file (e.g., DCD, XTC).
    trajrmsdstr : str
        Atom selection string for trajectory atoms used in RMSD calculation.
    trajalnstr : str or None, optional
        Atom selection string used to align the trajectory. If `None`, no
        explicit alignment selection is provided to MetricRmsd.
    refrmsdstr : str, optional
        Atom selection string for reference structure atoms
        (default: "protein and name CA").
    refalnstr : str, optional
        Atom selection string used to align the reference structure
        (default: "protein and name CA").
    centerstr : str, optional
        Atom selection string around which to center the wrapping of the
        trajectories (default: "protein and name CA").

    Output
    -------
    - CSV file: `../results/tables/{filename}_rmsd.csv`
     """
    csv_file = f"../results/tables/{filename}_rmsd.csv"
    
    if os.path.exists(csv_file):
        # If CSV exists, just load it to confirm readability and exit.
        pd.read_csv(csv_file)
    else:
        print('If not removed yet, run "filter_prod.py" to remove water before loading.')
        mol = Molecule(topology, validateElements=False)
        mol.read(trajectory)
        # mol.filter("not resname TIP3")  # run filter_prod.py before calling this function
        ref_mol = mol.copy(frames=[0])
        time_values = np.arange(0, mol.numFrames) / 10

        # Compute RMSD
        met = metricrmsd.MetricRmsd(
            ref_mol,
            trajrmsdstr=trajrmsdstr,
            trajalnstr=trajalnstr,
            refrmsdstr=refrmsdstr,
            refalnstr=refalnstr,
            centerstr=centerstr,
        )
        rmsd_values = met.project(mol)

        # Store data in CSV
        data = pd.DataFrame({
            'time': time_values,
            'rmsd': rmsd_values
        })
        data.to_csv(csv_file, index=False)

    return
