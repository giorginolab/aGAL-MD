# α-galactosidase MD preparation and analysis

The repo contains:

1. *glycosilation/*, contains the reglycosilated pdb structure and other glycan-related data.
2. *DGJ/*, contains DGJ chain A and B mol2 files used in mol_prep_fabry.py and other DGJ-related data.
3. *results/*, for each model and replica, contains csv tables and relative plots.
4. *functions/*, contains:
    - mol_prep_fabry.ipynb to prepare the system for equilibration;
    - check_end.py to quickly check in multiple sistems the equilibration/production termination;
    - production_prep.py to generate the production folder; 
    - evaluation.ipynb a notebook to compute and visualise RMSD and RMSF of different contitions (relies on rmsd.py and rmsf.py);
    - sbatch_acemd single file submission (either equilibration or production)
    - sbatch_many allow submission of consecutive files for long production (or equilibration).