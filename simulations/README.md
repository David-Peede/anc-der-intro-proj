# Simulations

This directory contains all the code to replicate the simulation based analyses from Peede et al. 202X, which I will now outline.

## Generate `demes` graphs as `.yaml` files.

To generate all `demes` graphs for all simulations run `create_yamls.sh` which will subsequently generate all `.yaml` files to run all analyses.

```bash
# Create all demes graphs.
./create_yamls.sh
```

## Run all simulations.

Now that all of the `demes` graphs are created, it is now time to generate the simulated data. In total this amounts to 17,200 simulation runs, which I broke up into multiple parts for convenience.

```bash
# Generate simulated data under the IUA model.
./distribute_simulations_part1.sh
./distribute_simulations_part2.sh
# Generate simulated data under the modfied Ragsdale and Gravel 2019 model.
./distribute_ragsdale_simulations_part1.sh
./distribute_ragsdale_simulations_part2.sh
# Generate simulated data under the original Ragsdale and Gravel 2019 model.
./distribute_ragsdale_simulations_org.sh
# Generate simulated data under the non-IUA models.
./distribute_non_iua_simulations_part1.sh
./distribute_non_iua_simulations_part2.sh
./distribute_non_iua_simulations_part3.sh
./distribute_non_iua_simulations_part4.sh
./distribute_non_iua_simulations_part5.sh
```

## Calculate observed site patterns and introgression statistics.

Next I calculated all observed site patterns and introgression statistics from each simulation. Again given the total number of simulations I broke this up into multiple parts for convenience.

```bash
# Calculate results for the IUA simulations.
./distribute_introgression_calculations_n1.sh
./distribute_introgression_calculations_n100_part1.sh
./distribute_introgression_calculations_n100_part2.sh
# Calculate results for the modfied Ragsdale and Gravel 2019 simulations.
./distribute_ragsdale_introgression_calculations_n1.sh
./distribute_ragsdale_introgression_calculations_n100_part1.sh
./distribute_ragsdale_introgression_calculations_n100_part2.sh
# Calculate results for the original Ragsdale and Gravel 2019 simulations.
./distribute_ragsdale_introgression_calculations_org.sh
# Calculate results for the non-IUA simulations.
./distribute_non_iua_introgression_calculations_n1.sh
./distribute_non_iua_introgression_calculations_n100_part1.sh
./distribute_non_iua_introgression_calculations_n100_part2.sh
./distribute_non_iua_introgression_calculations_n100_part3.sh
./distribute_non_iua_introgression_calculations_n100_part4.sh
./distribute_non_iua_introgression_calculations_n100_part5.sh
```

## Perform bootstrapping.

Lastly, for all replicate simulations under the IUA and Ragsdale and Gravel models each simulated genome needs to be bootstrapped 1000 times. Again given the amount of bootstrapping needed to be performed I broke this up into multiple parts for convenience.

```bash
# Generate bootstrapped distrubtuions for the IUA simulations.
./distribute_bootstraps_n1_part1.sh
./distribute_bootstraps_n1_part2.sh
./distribute_bootstraps_n100_part1.sh
./distribute_bootstraps_n100_part2.sh
# Generate bootstrapped distrubtuions for the modfied Ragsdale and Gravel 2019 simulations.
./distribute_ragsdale_bootstraps_n1_part1.sh
./distribute_ragsdale_bootstraps_n1_part2.sh
./distribute_ragsdale_bootstraps_n100_part1.sh
./distribute_ragsdale_bootstraps_n100_part2.sh
# Generate bootstrapped distrubtuions for the original Ragsdale and Gravel 2019 simulations.
./distribute_ragsdale_bootstraps_org.sh
```

## Analysis

A walkthrough of my main analysis for the IUA and Ragsdale and Gravel models can be viewed in the `main_simulation_analyses_v04.ipynb` notebook. A walkthrough of my supplemental analysis for the difference of site pattern differences on non-IUA models can be viewed in the `supplement_simulation_analyses_v04.ipynb` notebook.

## Notes

All shell scripts contain for loops to submit simulations in parallel using the `SLURM` scheduler on `OSCAR` at Brown University. You will need to modify the `SLURM` headers to run on your own HPC. I would recommend not changing the resource allocation for each job, as they have been optimized over numerous trials and errors.