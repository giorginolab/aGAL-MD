#!/usr/bin/env python3
"""Check whether ACEMD simulations completed successfully for each replica.

Run this script from the parent folder (e.g., within `3GXT_reglyco/`):

    python ../functions/check_end.py equilibration
    python ../functions/check_end.py production
"""
import sys
import os
import re

if len(sys.argv) < 2:
    print("Usage: check_end.py equilibration|production")
    sys.exit(1)
line_param = sys.argv[1]  # specify whether to check equilibration or production

replica_folders = sorted([d for d in os.listdir('.') if os.path.isdir(d)])
if not replica_folders:
    print("No replica folders found in current directory.")
    sys.exit(1)

pattern = re.compile(r"slurm-(\d+)\.out$")  # retrieve the submission number
for replica in replica_folders:
    target_folder = os.path.join(replica, line_param)
    if not os.path.isdir(target_folder):
        continue

    out_files = [f for f in os.listdir(target_folder) if f.endswith('.out')]

    max_num = -1
    max_file = None
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
