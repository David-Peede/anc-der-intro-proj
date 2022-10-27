### Dependencies ###
import msprime
import numpy as np
import sys
### sys.argv[1] = the neanderthal admixture proportion ###
### sys.argv[2] = the denisovan admixture proportion ###
### sys.argv[3] = replicate ID ###


# Define a two pulse model of introgression.
def two_pulse_human_model(admix_prop_nea, admix_prop_den):
    # Intialize demographic model.
    two_pulse_model = msprime.Demography()
    # We assume constant and equal effective population sizes for
    # all lineages.
    two_pulse_model.add_population(name='AFR', initial_size=10_000)
    two_pulse_model.add_population(name='EUR', initial_size=10_000)
    two_pulse_model.add_population(name='NEA', initial_size=10_000)
    two_pulse_model.add_population(name='DEN', initial_size=10_000)
    two_pulse_model.add_population(name='ARC', initial_size=10_000)
    two_pulse_model.add_population(name='AMH', initial_size=10_000)
    two_pulse_model.add_population(name='HUM', initial_size=10_000)
    # Introgression from the Denisovan to the Eurasian lineage
    # occuring 1,500 generations ago with a probability of admix_prop_den.
    two_pulse_model.add_mass_migration(
        time=1_500, source='EUR', dest='DEN', proportion=admix_prop_den,
    )
    # Introgression from the Neanderthal to the Eurasian lineage
    # occuring 2,200 generations ago with a probability of admix_prop_nea.
    two_pulse_model.add_mass_migration(
        time=2_200, source='EUR', dest='NEA', proportion=admix_prop_nea,
    )
    # The African and Eurasian lineages merge into the anatomically
    # modern human lineage 4,400 generations ago.
    two_pulse_model.add_population_split(
        time=4_400, derived=['AFR', 'EUR'], ancestral='AMH',
    )
    # The Denisovan and Neanderthal lineages merge into the archaic
    # lineage 17,080 generations ago.
    two_pulse_model.add_population_split(
        time=17_080, derived=['NEA', 'DEN'], ancestral='ARC',
    )
    # The anatomically modern human and archaic lineages merge
    # into the ancestral human lineage 26,320 generations ago.
    two_pulse_model.add_population_split(
        time=26_320, derived=['AMH', 'ARC'], ancestral='HUM',
    )
    return two_pulse_model

# Define simulator function.
def run_100mb_two_pulse_sim(admix_prop_nea, admix_prop_den, seed):
    # Simulate a 100 Mb tree-sequence.
    ts = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(1, ploidy=1, population='AFR'),
            msprime.SampleSet(1, ploidy=1, population='EUR'),
            msprime.SampleSet(1, ploidy=1, population='NEA'),
            ],
        demography=two_pulse_human_model(admix_prop_nea, admix_prop_den),
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
        './two_pulse/{0}_{1}/geno_mats/rep_id_{2}_geno_mat.csv.gz'.format(admix_prop_nea, admix_prop_den, seed),
        genotype_matrix,
        fmt='%d',
        delimiter=',',
        )
    # Extract the variable positions.
    variable_positions = mts.tables.sites.position
    # Save the variable positions.
    np.savetxt(
        './two_pulse/{0}_{1}/var_pos/rep_id_{2}_var_pos.csv.gz'.format(admix_prop_nea, admix_prop_den, seed),
        [variable_positions],
        fmt='%1.15f',
        delimiter=',',
        newline='\n',
    )
    # Save the tree-sequence just for good measures.
    mts.dump('./two_pulse/{0}_{1}/mut_tree_seq/rep_id_{2}_mut_tree_seq.ts'.format(admix_prop_nea, admix_prop_den, seed))
    return


# Parse comand-line arguments.
f_nea         = float(sys.argv[1])
f_den         = float(sys.argv[2])
rep_id    = int(sys.argv[3])

# Run the simulation!
run_100mb_two_pulse_sim(admix_prop_nea=f_nea, admix_prop_den=f_den, seed=rep_id)