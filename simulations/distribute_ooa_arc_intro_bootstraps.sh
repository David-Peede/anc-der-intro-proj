#!/bin/bash

for rep_id in {1..100}; do
sbatch -J bs_ooa_arc_intro_rep_id_${rep_id} -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/bs_ooa_arc_intro_rep_id_${rep_id}-%A.out -e logs/bs_ooa_arc_intro_rep_id_${rep_id}-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python generate_ooa_arc_intro_bootstrapped_data.py ${rep_id}"
done