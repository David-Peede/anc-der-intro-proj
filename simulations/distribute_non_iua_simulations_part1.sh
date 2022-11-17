#!/bin/bash

# I run these simulations in two parts since Brown will only allow 1.2K jobs to sit in the queue.
for rep_id in {1..100}; do for f_nea in 0.0 0.005 0.01 0.015 0.02; do
sbatch -J multi_${f_nea}_0.0_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/multi_${f_nea}_0.0_rep_id_${rep_id}-%A.out -e logs/multi_${f_nea}_0.0_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_non_iua_simulated_data.py ${f_nea} 0.0 multi ${rep_id} 1"
sbatch -J basal_${f_nea}_0.0_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/basal_${f_nea}_0.0_rep_id_${rep_id}-%A.out -e logs/basal_${f_nea}_0.0_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_non_iua_simulated_data.py ${f_nea} 0.0 basal ${rep_id} 1"
done; done

# I run these simulations in two parts since Brown will only allow 1.2K jobs to sit in the queue.
for rep_id in {1..100}; do for f_nea in 0.0 0.005 0.01 0.015 0.02; do
sbatch -J multi_${f_nea}_0.0_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/multi_${f_nea}_0.0_rep_id_${rep_id}-%A.out -e logs/multi_${f_nea}_0.0_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_non_iua_simulated_data.py ${f_nea} 0.0 multi ${rep_id} 100"
sbatch -J basal_${f_nea}_0.0_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/basal_${f_nea}_0.0_rep_id_${rep_id}-%A.out -e logs/basal_${f_nea}_0.0_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_non_iua_simulated_data.py ${f_nea} 0.0 basal ${rep_id} 100"
done; done