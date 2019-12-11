#!/bin/bash
#SBATCH --nodes=1 --ntasks-per-node=6
#SBATCH -p batch
source /home/oakley/.bashrc

infile=$1
treefile=$2
outfile=$3

hyphy meme --alignment $infile --tree $treefile --output $outfile
