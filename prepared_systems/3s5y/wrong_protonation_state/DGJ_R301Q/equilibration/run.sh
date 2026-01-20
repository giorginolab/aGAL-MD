#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate htmd
acemd >log.txt 2>&1