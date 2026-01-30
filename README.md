# α-galactosidase MD preparation and analysis

Development of molecular dynamics-based features from all-atom simulations of αGAL mutants to assess pharmacoperone responsiveness

|    |    |
| ----- | ------- |
| Preprint and code: | [![DOI](https://zenodo.org/badge/1067749794.svg)](https://doi.org/10.5281/zenodo.17651899) |
| Data: | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17552241.svg)](https://doi.org/10.5281/zenodo.17552241) |


## Summary
The repo is organized as follows:

- [*DGJ*](DGJ), contains mol2 files used to build the correct protein structure and corresponding forcefield parameters (charmm36).
- [*functions*](functions), contains notebooks and scripts used to obtain the data described in the report.
- [*glycosylation*](glycosylation), contains the [reglycosylated pdb structure](/glycosylation/3s5y_reglyco.pdb)* and information about the glycans.
- [*paper*](paper), this preprint: I. Cazzaniga, T. Giorgino, *Development of molecular dynamics-based features from all-atom simulations of αGAL mutants to assess pharmacoperone responsiveness*, 2025 ([preprint](report/aGAL_report.pdf)).
- [*prepared_systems*](prepared_systems), contains the six systems described in the report (apo and holo, wt, N215S, R301Q) before equilibrating.
- [*results*](results), contains csv tables and plots obtained from the analysis of the trajectories.

*NOTE* the repo lacks of a *dist* folder, where simulations are stored locally. 
 
## Quickstart (protocol overview)
For a more thorought explaination, check the [Tutorial](Tutorial) section below.

1. Install HTMD, Moleculekit and ACEMD software in a conda environment, then from the repository root run the system-preparation notebook:
    ```bash
        conda activate ace_software #leave it always active
        jupyter notebook functions/1_mol_prep.ipynb
    ```
   Running this notebook generates one `build/` and one `equilibration/` directory for each system (apo/holo and wild-type/mutant combinations) under `prepared_system/<pdb>`.
   *NOTE:* The notebook automatically generates either the wild-type system or the specified mutants, according to the configuration defined within the notebook.

2. On your HPC system, run equilibration and production:
    - generate the replicas of each system with [`functions/make_replicas.py`](functions/make_replicas.py). *NOTE:* replicas are stored in `dist/`:

    ```bash
       cd prepared_systems/<pdb>
       python ../../functions/make_replicas.py <n> #n is the number of replicas
    ```
    - submit all the equilibrations with [`sbatch_acemd`](functions/sbatch_acemd).
    - generate the `production folder/` for each equilibrated system
    - submit all the productions with [`sbatch_many`](functions/sbatch_many).

*NOTE:* use [`check_end.py`](functions/check_end.py) to verify the end of equilibrations or productions:
  
   
```bash
    cd dist/<pdb>
    python ../../check_end.py production #it is required to specify which folder to check
```

3. (Optional) strip water and subsample trajectories:
  
    ```bash
        cd dist/<pdb>
        python ../../functions/filter_prod.py
    ```
4. perform the standardized analysis:
   
    ```bash
        cd functions/
        jupyter notebook functions/2_evaluation.ipynb  # compute RMSD/RMSF tables and plots

        cd  results/<pdb>/tables
        python ../../../functions/residence_time.py    # compute ligand residence times
    ```
   The resulting CSV tables and plots are stored at `results/<pdb>/tables` and `results/<pdb>/plots`.

*The glycans in this pdb have been processed as described at step 0.

Our trajectories obtained from the PDB 3GXT systems in [`prepared_systems`](prepared_systems) folder are available at doi:[10.5281/zenodo.17552241](https://zenodo.org/records/17552241): 1 ns per frame, 3 replicas, water stripped out.

## Tutorial
This package relies on HTMD, Moleculekit and ACEMD software, which can be installed as:
```bash
    conda create -n ace_software #create a new conda environment (recommended but not necessary)
    conda activate ace_software #activate the new conda environment if you decided to create one
    conda install htmd acemd cuda-version=12 python=3.10 -c acellera -c conda-forge
```

### Step 0: reglycosylate the protein (optional)
The protein structure we used here is based on the [PDB 3S5Y](rcsb.org/structure/3S5Y/) model.
The original PDB underwent a reglycosylation step via [glycoshape.org](https://glycoshape.org/reglyco), to have the same glycan structure in each glycosylated sites, in particular (glycan ID G00026MO):
![Alt Text](glycosylation/glycan.svg)

### Step 1: system building and preparation
It can be easily done by running the [1_mol_prep.ipynb](functions/1_mol_prep.ipynb) notebook, which includes:

- mutation and glycosylation handling
- ligand insertion 
- system segmentation and preparation
- equilibration folder creation

Although the example workflow focuses on migalastat, other ligands or competitive inhibitors can be, in principle, modeled by providing the correct coordinates to fit in the dimer and by supplying compatible force-field parameters (e.g. via CGenFF), without changing the rest of the protocol.

The notebook generates a folder for each wild-type/mutant apo/holo structure, each generated folder contains a `build` folder with the intermediate steps of the system preparation and a `equilibration` folder. 
The folders generated are currently stored in [`prepared_systems`](prepared_system/3s5y).

If multiple replicas of the system are to be run, use [`make_replicas.py`](functions/make_replicas.py) to generate the desired number of replicas. 
   ```bash
        cd prepared_systems/3s5y
        python ../../functions/make_replicas.py <n> #n is the number of replicas desired, e.g. 3
   ```
 make_replicas.py currently generates a folder `dist` in which the replicas and the simulation data are stored. This folder and its content is not included in this public repository because of GitHub limitations of file size.

The local repo organization for the should be: 

```
    github_directory/
    └── prepared_systems/ 
    |   └── <pdb number>/
    |       └── <structure>/ (e.g. apo)
    |           ├── build/
    |           └── equilibration/
    └── dist/
    |   └── <pdb number>/
    |       └── <structure_replica>/ (e.g. apo_1, apo_2)
    |           ├── build/
    |           ├── equilibration/
    |           └── production/
    all the other folders
```

The following steps of this tutorial are carried within `dist` folder.

### Step 2: equilibration run
The equilibration must be run on HPC.
For each structure and replica contained in `dist`:

```bash
    conda activate ace_software
    
    cd dist/<pdb_number>/<structure_replica>/equilibration/
    sbatch ../../../../functions/sbatch_acemd #single submission, from the equilibration folder

    cd dist/<pdb_number>
    for dir in *; do [ -d "$dir/equilibration" ] || { echo "Skipping $dir: no equilibration folder"; continue; }; echo "Entering $dir/equilibration"; (cd "$dir/equilibration" && sbatch ../../../functions/sbatch_acemd); done #multiple submissions, from the <pdb_number>/ folder
```
**NOTE** the single equilibration submission must be launched from within the specific `equilibration` folder, the multiple submission must be launched from the `<pdb_number>` folder.

It is possible to check if the equilibrations are ended by running [check_end.py](functions/check_end.py), called from `<pdb_number>` folder:

```bash
    conda activate ace_software

    cd dist/<pdb_number>    
    python ../../functions/check_end.py equilibration
```
**NOTE** this file looks for the slurm file and checks if it ended, more precise evaluation must be carried.

### Step 3: production preparation and run
Once the equilibration is completed it is possible to generate the required files for the production run.

From within the specific `<structure_replica>` folder do:

```bash
    conda activate ace_software
    
    cd dist/<pdb_number>/<structure_replica>
    python ../../functions/production_prep.py
```
If running multiple systems in parallel, from the `<pdb_number>` folder, do:

```bash
    conda activate ace_software
    
    cd dist/<pdb_number>
    for dir in */; do [ -d "$dir" ] || continue; echo "Entering $dir"; (cd "$dir" && python ../../functions/production_prep.py); done
```
At this point, in both cases, it is possible to run the production. 
In some cases, the production run may require multiple job submissions due to limited walltime. To avoid manual resubmission, you can use [sbatch_many](functions/sbatch_many)
 to generate a chained sequence of SLURM job submissions.

```bash
    conda activate ace_software 
    
    #single structure/replica submission
    cd dist/<pdb_number>/<structure_replica>
    sbatch sbatch_many <n> 

    #multiple structure/replica submission
    cd dist/<pdb_number>
    for dir in *; do [ -d "$dir/production" ] || ( echo "Skipping $dir: no production folder"; continue; ); echo "Entering $dir/production"; (cd "$dir/production" && sbatch ../../../functions/sbatch_many <n>); done #multiple submissions in parallel
```
With *<n>* being the number of jobs to be generated for each structure_replica (for example, in our case, with a 1 μs long simulation and 24 h walltime, we set n = 7).

If multiple systems are run in parallel, it is possible to check if the equilibrations are ended by running [check_end.py](functions/check_end.py) from within the parent folder generated by the notebook:

```bash
    cd dist/<pdb_number>
    
    python ../../functions/check_end.py production
```
**NOTE** this file looks for the slurm file with the highest number (the last run) and checks if it ended correctly, more precise evaluation must be carried.

### Step 4: 
We included a set of simple analysis based on RMSD and RMSF on protein, ligand and glycans. 
This involves:

1. An optional [filter_prod.py](functions/filter_prod.py) filtering step to remove water (TIP3, otherwise needs to be fixed) and reduce the size of both the topology (psf) and trajectory (xtc). 
    The new trajectory is generated by skipping every 10 frames.
        This can be run from within the `<pdb_number>` folder as:
```bash
    cd dist/<pdb_number>

    python ../../functions/filter_prod.py
```

2. A notebook for computing the RMSD and RMSF analysis, [2_evaluation.ipynb](functions/2_evaluation.ipynb), with some examples we included in our report. 
    Both rmsd and rmsf functions included in the notebook generate csv files and plots stored in the [results/](results/) folder.

3. Python script to compute the residence time of each DGJ in the corresponding monomer, called [residence_time.py](functions/residence_time.py). 
    This can be run from within the `functions` folder as:
```bash
    cd  results/<pdb>/tables

    python ../../../functions/residence_time.py 
```
## SWITCHING TO AMBER FF
Amber requires a different labelling of the glycans compared to CHARMM-36 ff patch system.

### Step 1: glycosylate the protein
Starting from the original pdb structure, we attach the correct glycan to the N- sites (N139, N192 and N215) using [GLYCAM-WEB](https://glycam.org/gp/)
The glycan string is: DManpa1-3[DManpa1-6]DManpb1-4DGlcpNAcb1-4DGlcpNAcb1-OH .
This pdb is stored as [3s5y_amber.pdb](glycosylation/3s5y_amber.pdb).

### Step 2: generate DGJ-specific forcefield
    Generate the DGJ-specific ff files using antechamber and parmchk2.
```bash
    conda activate ace_software

    cd DGJ/

    antechamber -i DGJ_A_cgeneff.mol2 -fi mol2 -o dgj_amber.prepi -fo prepi -c bcc -nc 1 -rn DGJ -at gaff2
    parmchk2 -i dgj_amber.prepi -f prepi -o dgj.frcmod 
```    
   The forcefield works for both DGJ_A and DGJ_B.

### Step 3: fix output pdb and build the system for equilibration
GLYCAM-WEB outputs a pdb where resids are reordered (chian A starts at 1 instead of 32, chain B at 391 instead of 32).
To fix this issue and allow correct recognition of the N-glycosilated sites, run:

```bash
    conda activate ace_software

    python functions/fix_glycam_pdb.py
```
After that it is possible to:
```bash
        conda activate ace_software 
        jupyter notebook functions/mol_prep_amber.ipynb
```

From this point on, the steps are identical to the ones explained under the Tutorial section.

*NOTA* al momento abbiamo provato a fare run senza glicani, per aggiungerli, lasciare la lista resnames_to_remove = [] vuota (seconda cella del notebook)

## Acknowledgements

The report was conducted as part of the PROPHECY-GlycoRare project, funded by Partenariato Esteso “Health Extended ALliance for Innovative Therapies, Advanced Lab-research, and Integrated Approaches of Precision Medicine -- HEAL ITALIA -- PE 00000019”, a valere sulle risorse del Piano Nazionale di Ripresa e Resilienza (PNRR) Missione 4 “Istruzione e Ricerca” – Componente 2 “Dalla Ricerca all'Impresa” -- Investimento 1.3, finanziato dall’Unione europea – NextGenerationEU -- a valere sull’Avviso pubblico del Ministero dell'Università e della Ricerca (MUR) n. 341 del 15.03.2022 (CUP Spoke leader: Università degli Studi di Milano-Bicocca – CUP H43C22000830006 – Spoke 5 “Next-Gen Therapeutics”).

    
