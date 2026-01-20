'''
Set up production simulations from an equilibration directory.

Before running, activate the environment:

    conda activate ace_software

    python ../../../functions/production_prep.py
'''

from acemd.protocols import setup_production
import os

def main():

    folder_work = os.getcwd()           # Use current working directory
    prod_run = '1 us'                   # Default production run length
    prod_temp = 300                     # Default temperature

    equil_dir = os.path.join(folder_work, "equilibration")
    prod_dir = os.path.join(folder_work, "production")

    setup_production(
        equildir=equil_dir,
        outdir=prod_dir,
        run=prod_run,
        temperature=prod_temp
    )

if __name__ == "__main__":
    main()


