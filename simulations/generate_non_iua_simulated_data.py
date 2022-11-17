### Dependencies ###
import demes
import msprime
import numpy as np
import sys
### sys.argv[1] = Neanderthal admixture proportion ###
### sys.argv[2] = Denisovan admixture proportion ###
### sys.argv[3] = non-iua model ###
### sys.argv[4] = replicate ID ###
### sys.argv[5] = sample size ###


# Define a function to run a 100 Mb simulation of a non-iua model.
def run_100mb_non_iua_sim(admix_prop_nea, admix_prop_den, model, seed, n):
    # Define a sample size dictionary.
    n_dicc = {
        1: {'AFR': 1, 'EUR': 1, 'ASN': 1, 'NEA': 1},
        100: {'AFR': 100, 'EUR': 100, 'ASN': 100, 'NEA': 1},
    }
    # Define a model dictionary.
    model_dicc = {
        'multi': './yamls/multi_',
        'basal': './yamls/basal_',
    }
    # Load the demes graph.
    non_iua_graph = demes.load(model_dicc[model]+'nea_{0}_den_{1}.yaml'.format(admix_prop_nea, admix_prop_den))
    # Simulate a 100 Mb tree-sequence.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(n_dicc[n]['AFR'], ploidy=1, population='AFR'),
            msprime.SampleSet(n_dicc[n]['EUR'], ploidy=1, population='EUR'),
            msprime.SampleSet(n_dicc[n]['ASN'], ploidy=1, population='ASN'),
            msprime.SampleSet(n_dicc[n]['NEA'], ploidy=1, population='NEA'),
            ],
        demography=msprime.Demography.from_demes(non_iua_graph),
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
        './non_iua_models/{0}/n_{1}/{2}_{3}/geno_mats/rep_id_{4}_geno_mat.csv.gz'.format(model, n, admix_prop_nea, admix_prop_den, seed),
        genotype_matrix,
        fmt='%d',
        delimiter=',',
        )
    # Extract the variable positions.
    variable_positions = mts.tables.sites.position
    # Save the variable positions.
    np.savetxt(
        './non_iua_models/{0}/n_{1}/{2}_{3}/var_pos/rep_id_{4}_var_pos.csv.gz'.format(model, n, admix_prop_nea, admix_prop_den, seed),
        [variable_positions],
        fmt='%1.15f',
        delimiter=',',
        newline='\n',
    )
    # Save the tree-sequence just for good measures.
    mts.dump('./non_iua_models/{0}/n_{1}/{2}_{3}/mut_tree_seq/rep_id_{4}_mut_tree_seq.ts'.format(model, n, admix_prop_nea, admix_prop_den, seed))
    return


# Parse comand-line arguments.
f_nea     = float(sys.argv[1])
f_den     = float(sys.argv[2])
non_iua   = str(sys.argv[3])
rep_id    = int(sys.argv[4])
samp_size = int(sys.argv[5])

# Run the simulation!
run_100mb_non_iua_sim(admix_prop_nea=f_nea, admix_prop_den=f_den, model=non_iua, seed=rep_id, n=samp_size)