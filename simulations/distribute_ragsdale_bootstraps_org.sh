#!/bin/bash

# N =1.
for rep_id in {1..100}; do
sbatch -J bs_ragsdale_org_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=2G --account=ccmb-condo -o logs/bs_ragsdale_org_rep_id_${rep_id}-%A.out -e logs/bs_ragsdale_org_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_ragsdale_bootstrapped_data_n1_org.py ${rep_id}"
done

# N =100.
for rep_id in {1..100}; do
sbatch -J bs_ragsdale_org_rep_id_${rep_id} -N 1 -n 10 -t 1:00:00 --mem=4G --account=ccmb-condo -o logs/bs_ragsdale_org_rep_id_${rep_id}-%A.out -e logs/bs_ragsdale_org_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_ragsdale_bootstrapped_data_n100_org.py ${rep_id}"
done