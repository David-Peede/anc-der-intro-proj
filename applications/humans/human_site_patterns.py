### Dependencies ###
import allel
import numpy as np
import numcodecs
import pandas as pd
import random
import sys
import time
import zarr


# Define a function to subset populations from the TGP genotype matrix.
def pop_subset(
    meta_info,
    geno_mat,
    pop,
):
    """
    ###########################################################################
    INPUT
        meta_info: Metadata dataframe.
        geno_mat: Original genotype matrix.
        pop: TGP population to subset.
    ---------------------------------------------------------------------------
    OUTPUT:
        alt_freq_array: Alternative allele frequency array for the population.
    ###########################################################################
    """
    # Identify individuals in the target population.
    vals = meta_info['pop'].isin([pop]).values
    # Subset the meta info for only the target population.
    subset = meta_info[vals]
    # Get the index values of the indivuals in the target population.
    idx = subset.index.tolist()
    # Subset the genotype array for the target population and calculate alternative
    # allele frequencies.
    alt_freq_array = geno_mat.take(idx, axis=1).count_alleles().to_frequencies()[:, 1]
    return alt_freq_array

# Define a function to create genome wide alternative frequencies and positions dictionaries.
def genome_diccs():
    """
    ###########################################################################
    INPUT: N/A
    ---------------------------------------------------------------------------
    OUTPUT
        pop_dicc: Dictionary of alternative allele frequencies per chromosome,
                  per population.
        pos_dicc: Dictionary of variable position arrays per chromosome.
    ###########################################################################
    """
    # Intialize the populations to track.
    pop_list = [
        'BEB', 'STU', 'ITU', 'PJL', 'GIH',
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
        'PEL', 'MXL', 'CLM', 'PUR',
        'YRI', 'ALT', 'ANC',
    ]
    # Intialize empty dictionaries.
    pop_dicc = {}
    # For every population....
    for pop in pop_list:
        # Intialize a nested dictionary per chromosome.
        pop_dicc[pop] = {
            1: {}, 2: {}, 3: {}, 4: {}, 5: {},
            6: {}, 7: {}, 8: {}, 9: {}, 10: {},
            11: {}, 12: {}, 13: {}, 14: {},
            15: {}, 16: {}, 17: {}, 18: {},
            19: {}, 20: {}, 21: {}, 22: {}, 'X': {},
        }
    # Initialize empty dictionary to store variable position arrays.
    pos_dicc = {}
    # Intialize chromosome list.
    chrom_list = [
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14,
        15, 16, 17, 18,
        19, 20, 21, 22, 'X',
    ]
    # Read in the metadata.
    panel_path = './zarr_arrays/tgp_chimp_altai_info.txt'
    panel = pd.read_csv(
        panel_path, sep='\t', usecols=['sample', 'pop', 'super_pop']
    )
    # For every chromosome...
    for chrom in chrom_list:
        # Read in the zarr array and convert to a genotype matrix.
        zarr_path = './zarr_arrays/tgp_chimp_altai_merged_filtered_biallelic_chr{0}.zarr'.format(chrom)
        callset = zarr.open_group(zarr_path, mode='r')
        chrom_geno_mat = allel.GenotypeArray(callset['{0}/calldata/GT'.format(chrom)])
        # Build the positions dictionary.
        pos_dicc[chrom] = allel.SortedIndex(callset['{0}/variants/POS'.format(chrom)])
        # For every population...
        for pop in pop_list:
            # Build the population's genome dictionary.
            pop_dicc[pop][chrom] = pop_subset(panel, chrom_geno_mat, pop)
    return pop_dicc, pos_dicc

# Define a function to calculate site pattern counts.
def human_freq_site_patterns(
    P1,
    P2,
    P3,
    P4,
):
    """
    ###########################################################################
    INPUT
        P1: P1 whole-genome sequence.
        P2: P2 whole-genome sequence.
        P3: P3 whole-genome sequence.
        P4: P4 whole-genome sequence.
    ---------------------------------------------------------------------------
    OUTPUT: Genome counts of ABBA, BABA, BAAA, ABAA, and there HOM sites.
    ###########################################################################
    """
    # Initialize site pattern counts.
    ABBA = 0
    ABBA_HOM = 0
    BABA = 0
    BABA_HOM = 0
    BAAA = 0
    BAAA_HOM = 0
    ABAA = 0
    ABAA_HOM = 0
    # Loop through every variable site and calculate site pattern counts.
    for i in range(len(P4)):
        # When the ancestral allele is the reference allele.
        if P4[i] == 0:
            ABBA += ((1 - P1[i]) * P2[i] * P3[i] * (1 - P4[i]))
            ABBA_HOM += ((1 - P1[i]) * P3[i] * P3[i] * (1 - P4[i]))
            BABA += (P1[i] * (1 - P2[i]) * P3[i] * (1 - P4[i]))
            BABA_HOM += (P1[i] * (1 - P3[i]) * P3[i] * (1 - P4[i]))
            BAAA += (P1[i] * (1 - P2[i]) * (1 - P3[i]) * (1 - P4[i]))
            BAAA_HOM += (P1[i] * (1 - P3[i]) * (1 - P3[i]) * (1 - P4[i]))
            ABAA += ((1 - P1[i]) * P2[i] * (1 - P3[i]) * (1 - P4[i]))
            ABAA_HOM += ((1 - P1[i]) * P3[i] * (1 - P3[i]) * (1 - P4[i]))
        # When the ancestral allele is the alternative allele.
        elif P4[i] == 1:
            ABBA += (P1[i] * (1 - P2[i]) * (1 - P3[i]) * P4[i])
            ABBA_HOM += (P1[i] * (1 - P3[i]) * (1 - P3[i]) * P4[i])
            BABA += ((1 - P1[i]) * P2[i] * (1 - P3[i]) * P4[i])
            BABA_HOM += ((1 - P1[i]) * P3[i] * (1 - P3[i]) * P4[i])
            BAAA += ((1 - P1[i]) * P2[i] * P3[i] * P4[i])
            BAAA_HOM += ((1 - P1[i]) * P3[i] * P3[i] * P4[i])
            ABAA += (P1[i] * (1 - P2[i]) * P3[i] * P4[i])
            ABAA_HOM += (P1[i] * (1 - P3[i]) * P3[i] * P4[i])
    return ABBA, BABA, BAAA, ABAA, ABBA_HOM, BABA_HOM, BAAA_HOM, ABAA_HOM

# Create genome-wide alternative allele frequency and positions dictionaries.
tgp_genomes, _ = genome_diccs()
# Intialize the P2 populations.
p2_list = [
    'BEB', 'STU', 'ITU', 'PJL', 'GIH',
    'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
    'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
    'PEL', 'MXL', 'CLM', 'PUR',
]
# Intialize a dictionary to store bootstrapped values.
results_dicc = {}
# For all P2 populations.
for p2 in p2_list:
    # Fill the dictionary with empty site pattern lists.
    results_dicc[p2] = {
         'ABBA': 0, 'ABBA_HOM': 0,
         'BABA': 0, 'BABA_HOM': 0,
         'BAAA': 0, 'BAAA_HOM': 0,
         'ABAA': 0, 'ABAA_HOM': 0,
     }
# Intialize chromosome list.
chromosome_list = [
    1, 2, 3, 4, 5,
    6, 7, 8, 9, 10,
    11, 12, 13, 14,
    15, 16, 17, 18,
    19, 20, 21, 22, 'X',
]
# For ever P2 population....
for p2 in p2_list:
    # For all chromosomes.
    for chromosome in chromosome_list:
        # Calculate site pattern counts for the bootstrapped replicate and store the results.
        abba, baba, baaa, abaa, abba_hom, baba_hom, baaa_hom, abaa_hom = human_freq_site_patterns(
            P1=tgp_genomes[chromosome]['YRI'],
            P2=tgp_genomes[chromosome][p2],
            P3=tgp_genomes[chromosome]['ALT'],
            P4=tgp_genomes[chromosome]['ANC'],
        )
        results_dicc[p2]['ABBA'].append(abba)
        results_dicc[p2]['ABBA_HOM'].append(abba_hom)
        results_dicc[p2]['BABA'].append(baba)
        results_dicc[p2]['BABA_HOM'].append(baba_hom)
        results_dicc[p2]['BAAA'].append(baaa)
        results_dicc[p2]['BAAA_HOM'].append(baaa_hom)
        results_dicc[p2]['ABAA'].append(abaa)
        results_dicc[p2]['ABAA_HOM'].append(abaa_hom)
# Intialize a results file.
results_file = open('./tgp_site_patterns.csv', 'w')
# For all P2 populations...
for p2 in p2_pops:
    # Compile the results for that population.
    result_line = [
        p2,
        str(bs_dicc[p2]['ABBA']),
        str(bs_dicc[p2]['BABA']),
        str(bs_dicc[p2]['BAAA']),
        str(bs_dicc[p2]['ABAA']),
        str(bs_dicc[p2]['ABBA_HOM']),
        str(bs_dicc[p2]['BABA_HOM']),
        str(bs_dicc[p2]['BAAA_HOM']),
        str(bs_dicc[p2]['ABAA_HOM']),
    ]
    # Write the results line to the results file.
    results_file.write(','.join(result_line)+'\n')
# Close the results file.
results_file.close()