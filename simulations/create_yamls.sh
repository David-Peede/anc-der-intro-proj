#!/bin/bash

# Generate IUA yamls.
sbatch -J iua_yamls -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/iua_yamls-%A.out -e logs/iua_yamls-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_iua_yamls.py"

# Generate Ragsdale yamls.
sbatch -J ragsdale_yamls -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/ragsdale_yamls-%A.out -e logs/ragsdale_yamls-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_ragsdale_yamls.py"

# Generate Non-IUA yamls.
sbatch -J non_iua_yamls -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/non_iua_yamls-%A.out -e logs/non_iua_yamls-%A.err --mail-type=ALL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_non_iua_yamls.py"