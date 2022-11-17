#!/bin/bash

# Calculate the observed values.
for f in 0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5; do
sbatch -J obs_vals_${f} -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/obs_vals_${f}-%A.out -e logs/obs_vals_${f}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python calculate_introgression_statistics_n1.py ${f}"
done