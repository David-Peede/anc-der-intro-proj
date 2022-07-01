### Dependencies ###
import numpy as np
import pandas as pd
import random
import sys


# Define site pattern function.
def butterfly_site_patterns(
    p1,
    p2,
    p3,
):
    """
    ###########################################################################
    INPUT
        p1: Derived alle frequency array for the P1 lineage.
        p2: Derived alle frequency array for the P2 lineage.
        p3: Derived alle frequency array for the P3 lineage.
    ---------------------------------------------------------------------------
    OUTPUT: Genome counts of ABBA, BABA, BBAA, BAAA, ABAA, and AABA sites.
    ###########################################################################
    """
    # Calculate site pattern counts.
    abba = ((1 - p1) * (p2) * (p3)).sum()
    baba = ((p1) * (1 - p2) * (p3)).sum()
    baaa = ((p1) * (1 - p2) * (1 - p3)).sum()
    abaa = ((1 - p1) * (p2) * (1 - p3)).sum()
    return abba, baba, baaa, abaa

# Define bootstrapping function.
def butterfly_bootstrap(
    genome_df,
    p1_id,
    p2_id,
    p3_id,
    n_replicates,
):
    """
    ###########################################################################
    INPUT
        genome_df: Whole-genome dataframe for butterflies.
        p1_id: ID for the P1 population.
        p2_id: ID for the P2 population.
        p3_id: ID for the P3 population.
        n_replicates: Number of bootstrap replicates.
    ---------------------------------------------------------------------------
    OUTPUT: List of site pattern counts from bootstrapped replicates.
    ###########################################################################
    """
    # Intialize arrays to store bootstrapped values.
    abba_array = np.array([])
    baba_array = np.array([])
    baaa_array = np.array([])
    abaa_array = np.array([])
    abba_hom_array = np.array([])
    baba_hom_array = np.array([])
    baaa_hom_array = np.array([])
    abaa_hom_array = np.array([])
    # Intialize a dictionary of chromosome lengths from table S4.6.1 (doi:10.1038/nature11041).
    chromosome_dicc = {
        'chr1': 15_755_848, 'chr2': 3_590_614, 'chr3': 8_943_778,
        'chr4': 6_670_749, 'chr5': 8_000_463, 'chr6': 13_155_959,
        'chr7': 11_792_147, 'chr8': 6_828_437, 'chr9': 8_311_378,
        'chr10': 17_492_759, 'chr11': 11_233_900, 'chr12': 15_835_627,
        'chr13': 13_526_828, 'chr14': 6_527_076, 'chr15': 7_742_867,
        'chr16': 9_269_283, 'chr17': 13_935_249, 'chr18': 15_453_272,
        'chr19': 14_778_520, 'chr20': 5_748_428, 'chrZ': 4_088_377,
    }
    # For each bootstrap replicate.
    for i in range(n_replicates):
        # Intialize empty arrays to store bootstrapped genomes.
        bootstrapped_p1 = np.array([])
        bootstrapped_p2 = np.array([])
        bootstrapped_p3 = np.array([])
        # For each bootstrap window (219 x 1 Mb).
        for j in range(219):
            # Randomly select a chromosome.
            chromosome = random.choice(list(chromosome_dicc.keys()))
            # Randomly generate a start and end position.
            start = np.random.randint((chromosome_dicc[chromosome] - 999_999))
            end = start + 1_000_000
            # Subset the whole-genome dataframe to only include the randomly selected chromosome.
            chromosome_df = genome_df[genome_df['chromosome'] == chromosome]
            # Extract the variable positions for that chromosome.
            variable_positions = chromosome_df['chrom_pos'].values
            # Identify the variants that fall within the 1 Mb window.
            variants = np.where(((start <= variable_positions) & (variable_positions <= end)))[0]
            # Extract the focal populations to subset.
            p1_chromosome = chromosome_df[p1_id].values
            p2_chromosome = chromosome_df[p2_id].values
            p3_chromosome = chromosome_df[p3_id].values
            # Subset the 1 Mb window from each focal population.
            p1_subset = p1_chromosome[variants]
            p2_subset = p2_chromosome[variants]
            p3_subset = p3_chromosome[variants]
            # Append 219 windows of 1 Mb to make bootstrapped genomes.
            bootstrapped_p1 = np.append(bootstrapped_p1, p1_subset)
            bootstrapped_p2 = np.append(bootstrapped_p2, p2_subset)
            bootstrapped_p3 = np.append(bootstrapped_p3, p3_subset)
        # Calculate site pattern counts for the bootstrapped replicate and store the results.
        abba, baba, baaa, abaa = butterfly_site_patterns(
            p1=bootstrapped_p1,
            p2=bootstrapped_p2,
            p3=bootstrapped_p3,
        )
        abba_array = np.append(abba_array, abba)
        baba_array = np.append(baba_array, baba)
        baaa_array = np.append(baaa_array, baaa)
        abaa_array = np.append(abaa_array, abaa)
        abba_hom, baba_hom, baaa_hom, abaa_hom = butterfly_site_patterns(
            p1=bootstrapped_p1,
            p2=bootstrapped_p3,
            p3=bootstrapped_p3,
        )
        abba_hom_array = np.append(abba_hom_array, abba_hom)
        baba_hom_array = np.append(baba_hom_array, baba_hom)
        baaa_hom_array = np.append(baaa_hom_array, baaa_hom)
        abaa_hom_array = np.append(abaa_hom_array, abaa_hom)
        print('completed bootstrapped replicate {0}'.format(i))
    return abba_array, baba_array, baaa_array, abaa_array, abba_hom_array, baba_hom_array, baaa_hom_array, abaa_hom_array


# Parse comand-line arguments.
p1_pop = str(sys.argv[1])
p2_pop = str(sys.argv[2])
p3_pop = str(sys.argv[3])

# Load in the butterfly derived allele frequencies.
butterfly_df = pd.read_csv('./butterflies_wgs_filtered_der_freqs.csv')

# Preform bootstrapping.
bs_abba, bs_baba, bs_baaa, bs_abaa, bs_abba_hom, bs_baba_hom, bs_baaa_hom, bs_abaa_hom = butterfly_bootstrap(
    genome_df=butterfly_df,
    p1_id=p1_pop,
    p2_id=p2_pop,
    p3_id=p3_pop,
    n_replicates=1000,
)

# Save each bootstrapped distribution as a gzipped csv.
results_prefix = './bootstraps/{0}_{1}_{2}/'.format(p1_pop, p2_pop, p3_pop)

np.savetxt(
    results_prefix+'abba.csv.gz',
    [bs_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baba.csv.gz',
    [bs_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baaa.csv.gz',
    [bs_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'abaa.csv.gz',
    [bs_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'abba_hom.csv.gz',
    [bs_abba_hom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baba_hom.csv.gz',
    [bs_baba_hom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baaa_hom.csv.gz',
    [bs_baaa_hom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'abaa_hom.csv.gz',
    [bs_abaa_hom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
