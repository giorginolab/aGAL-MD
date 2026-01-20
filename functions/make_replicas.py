#copyes the folders stored at prepared_systems/<molecule>/, and creates n replicas as requested
#Usage: run from <molecule>/ as: python make_replicas.py <n_replicas>

from pathlib import Path
import shutil
import sys


if len(sys.argv) != 2:
    print("Usage: make_replicas.py <n_replicas>")
    sys.exit(1)

try:
    n_replicas = int(sys.argv[1])
    if n_replicas < 1:
        raise ValueError
except ValueError:
    print("<n_replicas> must be a positive integer")
    sys.exit(1)

prepared_dir = Path.cwd()              # prepared_systems/<system>
system_name = prepared_dir.name
project_root = prepared_dir.parents[1]

dist_root = project_root / "dist" / system_name
dist_root.mkdir(parents=True, exist_ok=True)

source_dirs = sorted(
    d for d in prepared_dir.iterdir()
    if d.is_dir()
)

if not source_dirs:
    print("No subdirectories found in current directory.")
    sys.exit(1)

print(f"Found {len(source_dirs)} folders to replicate")
print(f"Creating {n_replicas} replicas per folder")
print(f"Destination: {dist_root}")

for src_dir in source_dirs:
    for i in range(1, n_replicas + 1):
        replica_dir = dist_root / f"{src_dir.name}_{i}"

        if replica_dir.exists():
            print(f"[SKIP] {replica_dir} already exists")
            continue

        print(f"[CREATE] {replica_dir}")
        shutil.copytree(src_dir, replica_dir)




