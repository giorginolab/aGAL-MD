#!/usr/bin/env python3
#outside the protein folder (i.e. within 3GXT_reglyco/)
#python ../functions/check_end.py equilibration/production
import sys
import os
import re

# Check for correct usage
if len(sys.argv) < 2:
    print("Usage: script.py equilibration|production")
    sys.exit(1)

# Get the folder type to look into: "equilibration" or "production"
line_param = sys.argv[1]

# List all folders in the current directory (assumed replicas)
replica_folders = sorted([d for d in os.listdir('.') if os.path.isdir(d)])

if not replica_folders:
    print("No replica folders found in current directory.")
    sys.exit(1)

# Compile regex to capture the number after 'slurm-' and before '.out'
pattern = re.compile(r"slurm-(\d+)\.out$")

for replica in replica_folders:
    target_folder = os.path.join(replica, line_param)
    if not os.path.isdir(target_folder):
        continue

    out_files = [f for f in os.listdir(target_folder) if f.endswith('.out')]

    max_num = -1
    max_file = None
    # Find the slurm-*.out file with the highest number
    for f in out_files:
        match = pattern.match(f)
        if not match:
            continue
        num = int(match.group(1))
        if num > max_num:
            max_num = num
            max_file = f

    if max_file is None:
        continue

    file_path = os.path.join(target_folder, max_file)
    found = False

    # Open the file and search for the string "Simulation completed!"
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if "Simulation completed!" in line:
                    found = True
                    break
    except Exception as e:
        print(f"Could not read file {file_path}: {e}")
        continue

    status = "FOUND" if found else "NOT FOUND"
    print(f"[{replica}] {line_param}/{max_file} : Simulation completed! {status}")
