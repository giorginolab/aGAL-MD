#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate ace_software
acemd >log.txt 2>&1