### Dependencies ###
import demes
import msprime
import numpy as np
import sys
### sys.argv[1] = replicate ID ###


# Define a function to run a 100 Mb simulation of the Ragsdale and Gravel 2019 model.
def run_100mb_ooa_arc_intro_sim(seed):
    # Load the demes graph.
    ooa_arc_intro_graph = demes.load('./HomSap__OutOfAfricaArchaicAdmixture_5R19.yaml')
    # Simulate a 100 Mb tree-sequence.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(1, ploidy=1, population='YRI'),
            msprime.SampleSet(1, ploidy=1, population='CEU'),
            msprime.SampleSet(1, ploidy=1, population='CHB'),
            msprime.SampleSet(1, ploidy=1, population='Neanderthal'),
            ],
        demography=msprime.Demography.from_demes(ooa_arc_intro_graph),
        sequence_length=100_000_000, recombination_rate=10**-8,
        random_seed=seed,
    )
    # Overlay mutations on the tree-sequence.
    mts = msprime.sim_mutations(
        tree_sequence=ts, rate=1.5 * 10**-8,
        model='jc69', random_seed=seed,
        discrete_genome=False,
    )
    # Extract the genotype matrix.
    genotype_matrix = mts.genotype_matrix()
    # Save the genotype matrix.
    np.savetxt(
        './sim_outputs/ooa_arc_intro/geno_mats/rep_id_{0}_geno_mat.csv.gz'.format(seed),
        genotype_matrix,
        fmt='%d',
        delimiter=',',
        )
    # Extract the variable positions.
    variable_positions = mts.tables.sites.position
    # Save the variable positions.
    np.savetxt(
        './sim_outputs/ooa_arc_intro/var_pos/rep_id_{0}_var_pos.csv.gz'.format(seed),
        [variable_positions],
        fmt='%1.15f',
        delimiter=',',
        newline='\n',
    )
    # Save the tree-sequence just for good measures.
    mts.dump('./sim_outputs/ooa_arc_intro/mut_tree_seq/rep_id_{0}_mut_tree_seq.ts'.format(seed))
    return


# Parse comand-line arguments.
rep_id = int(sys.argv[1])

# Run the simulation!
run_100mb_ooa_arc_intro_sim(seed=rep_id)