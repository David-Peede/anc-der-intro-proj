#!/bin/bash

# Calculate introgression metrics.
for f_nea in 0.0 0.005 0.01 0.015 0.02; do for f_den in 0.0 0.005 0.01 0.015 0.02; do
sbatch -J multi_obs_vals_${f_nea}_${f_den} -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/multi_obs_vals_${f_nea}_${f_den}-%A.out -e logs/multi_obs_vals_${f_nea}_${f_den}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python calculate_non_iua_introgression_statistics_n1.py ${f_nea} ${f_den} multi"
sbatch -J basal_obs_vals_${f_nea}_${f_den} -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/basal_obs_vals_${f_nea}_${f_den}-%A.out -e logs/basal_obs_vals_${f_nea}_${f_den}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python calculate_non_iua_introgression_statistics_n1.py ${f_nea} ${f_den} basal"
done; done
