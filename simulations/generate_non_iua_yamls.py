# Import packages.
import demes
import msprime


# Define a multi pulse model of introgression.
def multi_pulse_human_model(admix_prop_nea, admix_prop_den):
    # Intialize demographic model.
    multi_pulse_model = msprime.Demography()
    # We assume constant and equal effective population sizes for
    # all lineages.
    multi_pulse_model.add_population(name='AFR', initial_size=10_000)
    multi_pulse_model.add_population(name='EUR', initial_size=10_000)
    multi_pulse_model.add_population(name='ASN', initial_size=10_000)
    multi_pulse_model.add_population(name='EURS', initial_size=10_000)
    multi_pulse_model.add_population(name='NEA', initial_size=10_000)
    multi_pulse_model.add_population(name='DEN', initial_size=10_000)
    multi_pulse_model.add_population(name='ARC', initial_size=10_000)
    multi_pulse_model.add_population(name='AMH', initial_size=10_000)
    multi_pulse_model.add_population(name='HUM', initial_size=10_000)
    # Introgression from the Neanderthal to the Asian lineage
    # occuring 880 generations ago with a probability of admix_prop_nea.
    # (Modified from Jacobs et al 2019).
    multi_pulse_model.add_mass_migration(
        time=880, source='ASN', dest='NEA', proportion=admix_prop_nea,
    )
    # Introgression from the Denisovan to the Asian lineage
    # occuring 1,500 generations ago with a probability of admix_prop_den.
    # (Modified from Rogers et al 2015).
    multi_pulse_model.add_mass_migration(
        time=1_500, source='ASN', dest='DEN', proportion=admix_prop_den,
    )
    # The Asian and European lineages merge into the Eurasian
    # modern human lineage 1,600 generations ago.
    # (Modified from Rogers et al 2015).
    multi_pulse_model.add_population_split(
        time=1_600, derived=['ASN', 'EUR'], ancestral='EURS',
    )
    # Introgression from the Neanderthal to the Eurasian lineage
    # occuring 2,200 generations ago with a probability of admix_prop_nea.
    # (Modified from Sankararaman et al 2012).
    multi_pulse_model.add_mass_migration(
        time=2_200, source='EURS', dest='NEA', proportion=admix_prop_nea,
    )
    # The African and Eurasian lineages merge into the anatomically
    # modern human lineage 4,400 generations ago.
    # (Modified from Veeramah and Hammer 2014).
    multi_pulse_model.add_population_split(
        time=4_400, derived=['AFR', 'EURS'], ancestral='AMH',
    )
    # The Denisovan and Neanderthal lineages merge into the archaic
    # lineage 17,080 generations ago.
    # (Modified from Prufer et al 2014).
    multi_pulse_model.add_population_split(
        time=17_080, derived=['NEA', 'DEN'], ancestral='ARC',
    )
    # The anatomically modern human and archaic lineages merge
    # into the ancestral human lineage 26,320 generations ago.
    # (Modified from Prufer et al 2014).
    multi_pulse_model.add_population_split(
        time=26_320, derived=['AMH', 'ARC'], ancestral='HUM',
    )
    return multi_pulse_model

# Define a multi pulse model of introgression with gene flow from basal eurasians.
def basal_multi_pulse_human_model(admix_prop_nea, admix_prop_den):
    # Intialize demographic model.
    basal_multi_pulse_model = msprime.Demography()
    # We assume constant and equal effective population sizes for
    # all lineages.
    basal_multi_pulse_model.add_population(name='AFR', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='EUR', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='ASN', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='EURS', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='BASAL', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='ANC_EURS', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='NEA', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='DEN', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='ARC', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='AMH', initial_size=10_000)
    basal_multi_pulse_model.add_population(name='HUM', initial_size=10_000)
    # Introgression from the Neanderthal to the Asian lineage
    # occuring 880 generations ago with a probability of admix_prop_nea.
    # (Modified from Jacobs et al 2019).
    basal_multi_pulse_model.add_mass_migration(
        time=880, source='ASN', dest='NEA', proportion=admix_prop_nea,
    )
    # Dilution event from Basal Eurasians into Europeans
    # 1,000 generations ago with a probability of 0.2.
    # (Modified from Villenea and Schraiber 2019)
    basal_multi_pulse_model.add_mass_migration(
        time=1_100, source='EUR', dest='BASAL', proportion=0.2,
    )
    # Introgression from the Denisovan to the Asian lineage
    # occuring 1,500 generations ago with a probability of admix_prop_den.
    # (Modified from Rogers et al 2015).
    basal_multi_pulse_model.add_mass_migration(
        time=1_500, source='ASN', dest='DEN', proportion=admix_prop_den,
    )
    # The Asian and European lineages merge into the Eurasian
    # modern human lineage 1,600 generations ago.
    # (Modified from Rogers et al 2015).
    basal_multi_pulse_model.add_population_split(
        time=1_600, derived=['ASN', 'EUR'], ancestral='EURS',
    )
    # Introgression from the Neanderthal to the Eurasian lineage
    # occuring 2,200 generations ago with a probability of admix_prop_nea.
    # (Modified from Sankararaman et al 2012).
    basal_multi_pulse_model.add_mass_migration(
        time=2_200, source='EURS', dest='NEA', proportion=admix_prop_nea,
    )
    # The Basal and presenty day Eurasian lineages merge into the
    # Ancestral Eurasian lineage 3,000 generations ago.
    # (Modified from Villenea and Schraiber 2019).
    basal_multi_pulse_model.add_population_split(
        time=3_000, derived=['BASAL', 'EURS'], ancestral='ANC_EURS',
    )
    # The African and Ancestral Eurasian lineages merge into
    # the anatomically modern human lineage 4,400 generations ago.
    # (Modified from Veeramah and Hammer 2014).
    basal_multi_pulse_model.add_population_split(
        time=4_400, derived=['AFR', 'ANC_EURS'], ancestral='AMH',
    )
    # The Denisovan and Neanderthal lineages merge into the archaic
    # lineage 17,080 generations ago.
    # (Modified from Prufer et al 2014).
    basal_multi_pulse_model.add_population_split(
        time=17_080, derived=['NEA', 'DEN'], ancestral='ARC',
    )
    # The anatomically modern human and archaic lineages merge
    # into the ancestral human lineage 26,320 generations ago.
    # (Modified from Prufer et al 2014).
    basal_multi_pulse_model.add_population_split(
        time=26_320, derived=['AMH', 'ARC'], ancestral='HUM',
    )
    return basal_multi_pulse_model


# Define a list of admixture proportions.
nea_admix_props = [0.0, 0.005, 0.01, 0.015, 0.02]
den_admix_props = [0.0, 0.005, 0.01, 0.015, 0.02]

# For every Neanderthal admixture proportion...
for f_nea in nea_admix_props:
    # For every Denisovan admixture proportion...
    for f_den in den_admix_props:
        # Intialize a multi-pulse demographic model.
        multi_pulse_demo = multi_pulse_human_model(f_nea, f_den)
        # Convert the demographic model to a demes graph.
        multi_pulse_graph = msprime.Demography.to_demes(multi_pulse_demo)
        # Intialize a multi pulse with gene flow from basal eurasians demographic model.
        basal_demo = basal_multi_pulse_human_model(f_nea, f_den)
        # Convert the demographic model to a demes graph.
        basal_graph = msprime.Demography.to_demes(basal_demo)
        # Export the yamls.
        demes.dump(multi_pulse_graph, './yamls/multi_nea_{0}_den_{1}.yaml'.format(f_nea, f_den))
        demes.dump(basal_graph, './yamls/basal_nea_{0}_den_{1}.yaml'.format(f_nea, f_den))