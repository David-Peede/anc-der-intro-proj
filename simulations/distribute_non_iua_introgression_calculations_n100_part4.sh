#!/bin/bash

# I run these simulations in two parts since Brown will only allow 1.2K jobs to sit in the queue.
for rep_id in {1..100}; do for f_nea in 0.0 0.005 0.01 0.015 0.02; do
sbatch -J multi_calc_${f_nea}_0.015_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=3G --account=ccmb-condo -o logs/multi_calc_${f_nea}_0.015_rep_id_${rep_id}-%A.out -e logs/multi_calc_${f_nea}_0.015_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python calculate_non_iua_introgression_statistics_n100.py ${f_nea} 0.015 multi ${rep_id}"
sbatch -J basal_calc_${f_nea}_0.015_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=3G --account=ccmb-condo -o logs/basal_calc_${f_nea}_0.015_rep_id_${rep_id}-%A.out -e logs/basal_calc_${f_nea}_0.015_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python calculate_non_iua_introgression_statistics_n100.py ${f_nea} 0.015 basal ${rep_id}"
done; done