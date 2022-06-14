#!/bin/bash

# I run these simulations in two parts since Brown will only allow 1.2K jobs to sit in the queue.
for rep_id in {1..100}; do for f in 0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do
sbatch -J bs_${f}_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/bs_${f}_rep_id_${rep_id}-%A.out -e logs/bs_${f}_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_bootstrapped_data.py ${f} ${rep_id}"
done; done