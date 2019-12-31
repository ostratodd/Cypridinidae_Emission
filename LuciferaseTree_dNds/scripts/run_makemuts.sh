#!/bin/bash
sing_image="/home/pcw/ubuntu_bioperl_biopy2.img"
script_path="/home/oakley/projects/Cypridinidae_Emission/LuciferaseTree_dNds/scripts/makemutations.pl"
infile=$1
mutfile=$2
pheno=$3
singularity exec $sing_image perl $script_path $infile $mutfile $pheno
