# Import packages.
import demes
import msprime


# Define IUA model of introgression.
def iua_human_model(admix_prop):
    # Intialize demographic model.
    iua_model = msprime.Demography()
    # We assume constant and equal effective population sizes for
    # all lineages.
    iua_model.add_population(name='AFR', initial_size=10_000)
    iua_model.add_population(name='EUR', initial_size=10_000)
    iua_model.add_population(name='NEA', initial_size=10_000)
    iua_model.add_population(name='AMH', initial_size=10_000)
    iua_model.add_population(name='HUM', initial_size=10_000)
    # Introgression from the Neanderthal to the Eurasian lineage
    # occuring 1,600 generations ago with a probability of f.
    iua_model.add_mass_migration(
        time=1_600, source='EUR', dest='NEA', proportion=admix_prop,
    )
    # The African and Eurasian lineages merge into the anatomically
    # modern human lineage 4,000 generations ago.
    iua_model.add_population_split(
        time=4_000, derived=['AFR', 'EUR'], ancestral='AMH',
    )
    # The anatomically modern human and Neanderthal lineages merge
    # into the ancestral human lineage 16,000 generations ago.
    iua_model.add_population_split(
        time=16_000, derived=['AMH', 'NEA'], ancestral='HUM',
    )
    return iua_model



# Define a list of admixture proportions.
f_vals = [
    0.0, 0.01, 0.02,
    0.03, 0.04, 0.05,
    0.06, 0.07, 0.08,
    0.09, 0.1, 0.2,
    0.3, 0.4, 0.5,
]

# For every admixture proportion...
for f in f_vals:
    # Intialize an IUA demographic model.
    iua_demo = iua_human_model(f)
    # Convert the demographic model to a demes graph.
    iua_graph = msprime.Demography.to_demes(iua_demo)
    # Export the yamls.
    demes.dump(iua_graph, './yamls/iua_f_{0}.yaml'.format(f))