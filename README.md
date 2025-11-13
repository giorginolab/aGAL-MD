# α-galactosidase MD preparation and analysis

The repo contains:

- *DGJ*, contains mol2 files used to build the correct protein structure and corresponding forcefield parameters (charmm36).
- *functions*, contains notebooks and scripts used to obtain the data described in the report.
- *glycosylation*, contains the [reglycosylated pdb structure](/glycosylation/3gxt_reglyco.pdb)* and informations about the glycans.
- *prepared_systems*, contains the six systems described in the report (apo and holo, wt, N215S, R301Q) before equilibrating.
- *results*, contains csv tables and plots obtained from the analysis of the trajectories.
- *report*, a technical report in long-form.
 

*The glycans in this pdb have been processed as described at step 0.

Trajectories obtained from the molecules present in `prepared_systems` folder are available with DOI:[10.5281/zenodo.17552241](https://zenodo.org/records/17552241): 1/ns  fames, 3 replicas, water stripped out.

## Tutorial
This package relies on HTMD, Moleculekit and ACEMD software, which can be installed as:
```bash
    conda create -n ace_software #create a new conda environment (recommended but not necessary)
    conda activate ace_software #activate the new conda environment if you decided to create one
    conda install htmd acemd cuda-version=12 python=3.10 -c acellera -c conda-forge
```
### Step 0: reglycosylate the protein (NON mandatory)
The protein structure we used here is based on the [PDB 3GXT](rcsb.org/structure/3GXT/) model, which underwent a reglycosylation step via [glycoshape.org](https://glycoshape.org/reglyco) 
to have the same glycan structure in each glycosylated sites, in particular:
![Alt Text](glycosylation/glycan.svg)

The NOJ (ligand) will be addressed in step 1 of the tutorial.

### Step 1: system building and preparation
It can be easily done by running the [mol_prep_fabry.ipynb](functions/mol_prep_fabry.ipynb) which includes:

- mutation and glycosylation handling
- NOJ replacement with migalastat (DGJ)
- system segmentation and preparation
- equilibration folder creation

The notebook generates a folder for each mutant apo/holo structure, each folder contains a `build` folder with the intermediate steps of the system preparation and an `equilibration` folder. 

If multiple replicas of the system are to be run, please make copies of the folders at this point and store them in a parent folder following this pattern:

```
github_directory/
└── parent_folder/ (i.e. 3GXT_reglyco)
    └── structure_replica/ (i.e. apo_1)
        ├── build/
        └── equilibration/
```


### Step 2: equilibration run
The equilibration must be runned on HPC.
From within the `equilibration` folder:

```bash
    conda activate htmd 
    
    sbatch sbatch_acemd #single strucutre_replica

    for dir in *; do [ -d "$dir/equilibration" ] || { echo "Skipping $dir: no equilibration folder"; continue; }; echo "Entering $dir/equilibration"; (cd "$dir/equilibration" && sbatch ../../../functions/sbatch_acemd); done #multiple submissions
```
**NOTE** the single equilibration submission must be launched from within the specific `equilibration` folder, the multiple submission must be launched from the parernt folder.

If multiple systems are run in parallel, it is possible to check if the equilibrations are ended by running [check_end.py](functions/check_end.py), called from the parent folder:

```bash
    cd 3GXT_reglyco # change with the name of your parent folder
    
    python ../functions/check_end.py equilibration
```
**NOTE** this file looks for the slurm file and checks if it ended correctly, more precise evaluation must be carried.

### Step 3: production preparation and run
Once the equilibration is completed it is possible to generate the required files for the production run.

From within the specific structure_replica folder do:

```bash
    conda activate htmd
    
    python ../../functions/production_prep.py
```
If running multiple systems in parallel, from the parent folder, do:

```bash
    conda activate htmd
    
    for dir in */; do [ -d "$dir" ] || continue; echo "Entering $dir"; (cd "$dir" && python ../../functions/production_prep.py); done
```
At this point, in both cases, it is possible to run the production. 
In some cases the production could require multiple job submission if maximum runtime 24 hours, to avoid doing so manually we use [sbatch_many](functions/sbatch_many). 

```bash
    conda activate htmd 
    
    sbatch sbatch_many n  #single structure_replica
    for dir in *; do [ -d "$dir/production" ] || ( echo "Skipping $dir: no production folder"; continue; ); echo "Entering $dir/production"; (cd "$dir/production" && sbatch ../../../functions/sbatch_many n); done #multiple submissions in parallel
```
With *n* being the number of jobs to be generated for each structure_replica (for example, in our case a 1 μs long simulation required about 7 jobs).

If multiple systems are run in parallel, it is possible to check if the equilibrations are ended by running [check_end.py](functions/check_end.py) from within the parent folder generated by the notebook:

```bash
    cd 3GXT_reglyco #change with the name of your parent folder
    
    python ../functions/check_end.py production
```
**NOTE** this file looks for the slurm file with the highest number (the last run) and checks if it ended correctly, more precise evaluation must be carried.

### Step 4: 
We included a set of simple analysis based on RMSD and RMSF on protein, ligand and glycans. 
This involves:

1. A *non mandatory* [filter_prod.py](functions/filter_prod.py) filtering step to remove water (TIP3, otherwise needs to be fixed) and reduce the size of both the topology (psf) and trajectory (xtc). 
    The new trajectory is generated by skipping every 10 frames.
        This can be run from within the parent folder (i. e. 3GXT_reglyco) as:
```bash
        python residence_time.py
```

2. A notebook for computing the RMSD and RMSF analysis, [evaluation.ipynb](functions/evaluation.ipynb), with some examples we included in our report. 
    Both rmsd and rmsf functions included in the notebook generate csv files and plots stored in the [results/](results/) folder.

3. Python script to compute the residence time of each DGJ in the corresponding monomer, called [residence_time.py](functions/residence_time.py). 
    This can be run from within the `functions` folder as:
```bash
        python ../functions/filter_prod.py #run on GPUs in case of many structures
```


## Acknowledgements

The report was conduced as part of the  PROPHECY-GlycoRare project, funded by Partenariato Esteso “Health Extended ALliance for Innovative Therapies, Advanced Lab-research, and Integrated Approaches of Precision Medicine -- HEAL ITALIA -- PE 00000019”, a valere sulle risorse del Piano Nazionale di Ripresa e Resilienza (PNRR) Missione 4 “Istruzione e Ricerca” – Componente 2 “Dalla Ricerca all'Impresa” -- Investimento 1.3, finanziato dall’Unione europea – NextGenerationEU -- a valere sull’Avviso pubblico del Ministero dell'Università e della Ricerca (MUR) n. 341 del 15.03.2022 (CUP Spoke leader: Università degli Studi di Milano-Bicocca – CUP H43C22000830006 – Spoke 5 “Next-Gen Therapeutics”).

    
