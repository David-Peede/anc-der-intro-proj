#!/bin/bash

sbatch -J ooa_arc_no_intro_obs_vals -N 1 -n 1 -t 1:00:00 --mem=2G -o logs/ooa_arc_no_intro_obs_vals-%A.out -e logs/ooa_arc_no_intro_obs_vals-%A.err --mail-type=FAIL --mail-user=david_peede@brown.edu --wrap="module load anaconda/3-5.2.0; source /gpfs/runtime/opt/anaconda/3-5.2.0/etc/profile.d/conda.sh; conda activate anc_der_intro_env; python calculate_ooa_arc_no_intro_introgression_statistics.py"
