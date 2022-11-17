#!/bin/bash

# I run these simulations in two parts since Brown will only allow 1.2K jobs to sit in the queue.
for rep_id in {1..100}; do for f in 0.2 0.3 0.4 0.5; do
sbatch -J ragsdale_obs_vals_${f}_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=3G --account=ccmb-condo -o logs/ragsdale_obs_vals_${f}_rep_id_${rep_id}-%A.out -e logs/ragsdale_obs_vals_${f}_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python calculate_ragsdale_introgression_statistics_n100.py ${f} ${rep_id}"
done; done