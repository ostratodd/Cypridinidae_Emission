#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=5
#SBATCH -p short
source /home/oakley/.bashrc

iqtree -s /home/oakley/luciferasetree/sequences/clean_all_aa_ALN.fa -nt 5 -redo -pre /home/oakley/luciferasetree/results/phylogenies/all_aa

