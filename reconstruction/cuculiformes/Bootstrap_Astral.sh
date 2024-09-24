#!/bin/bash

#SBATCH --get-user-env
#SBATCH --mail-user=merondun@bio.lmu.de
#SBATCH --clusters=biohpc_gen
#SBATCH --partition=biohpc_gen_normal
#SBATCH --cpus-per-task=6
#SBATCH --time=12:00:00

java -jar ~/modules/ASTRAL/astral.5.7.7.jar -i hackett_sequenced_species.tre -o species_tree.tre --gene-only -r 2

