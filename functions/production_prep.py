from acemd.protocols import setup_production
import os

def main():
    """Set up production simulations from an equilibration directory."""
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

# Example usage from within a replica folder:
#   python ../../functions/production_prep.py
#
# To run across multiple replica folders from the parent directory:
#   for dir in *_2/; do
#       [ -d "$dir" ] || continue
#       echo "Entering $dir"
#       (cd "$dir" && python ../../functions/production_prep.py)
#   done
