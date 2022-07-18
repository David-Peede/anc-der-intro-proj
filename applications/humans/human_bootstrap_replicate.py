### Dependencies ###
import allel
import multiprocessing
import numpy as np
import numcodecs
import pandas as pd
import random
import sys
import time
import zarr


class Window:
    # Static variables, only one copy exists.
    pop_dicc = {}
    pos_dicc = {}

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

    # Define a function to create genome wide alternative frequencies and position dictionaries.
    def initialize_big_dictionaries():
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
        # For every population....
        for pop in pop_list:
            # Intialize a nested dictionary per chromosome.
            Window.pop_dicc[pop] = {
                1: {}, 2: {}, 3: {}, 4: {}, 5: {},
                6: {}, 7: {}, 8: {}, 9: {}, 10: {},
                11: {}, 12: {}, 13: {}, 14: {},
                15: {}, 16: {}, 17: {}, 18: {},
                19: {}, 20: {}, 21: {}, 22: {}, 'X': {},
            }
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
            Window.pos_dicc[chrom] = allel.SortedIndex(callset['{0}/variants/POS'.format(chrom)])
            # For every population...
            for pop in pop_list:
                # Build the population's genome dictionary.
                Window.pop_dicc[pop][chrom] = Window.pop_subset(panel, chrom_geno_mat, pop)

    # Define a function to work on an instance of the window class.
    def do_work(self):
        """
        ###########################################################################
        OUTPUT
         results_dicc: Dictionary of bootstrapped site patterns.
        ###########################################################################
        """
        # Track the window start time.
        win_start = time.time()
        # Intialize the P2 populations.
        p2_list = [
         'BEB', 'STU', 'ITU', 'PJL', 'GIH',
         'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
         'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
         'PEL', 'MXL', 'CLM', 'PUR',
        ]
        # Intialize a dictionary to store bootstrapped values.
        results_dicc = {}
        # For all P2 populations...
        for p2 in p2_list:
            # Fill the dictionary with empty site pattern lists.
            results_dicc[p2] = {
                 'ABBA': 0, 'ABBA_HOM': 0,
                 'BABA': 0, 'BABA_HOM': 0,
                 'BAAA': 0, 'BAAA_HOM': 0,
                 'ABAA': 0, 'ABAA_HOM': 0,
             }
        # Intialize a dictionary of chromosome lengths for hg19.
        chromosome_dicc = {
         1: 249250621, 2: 243199373, 3: 198022430,
         4: 191154276, 5: 180915260, 6: 171115067,
         7: 159138663, 8: 146364022, 9: 141213431,
         10: 135534747, 11: 135006516, 12: 133851895,
         13: 115169878, 14: 107349540, 15: 102531392,
         16: 90354753, 17: 81195210, 18: 78077248,
         19: 59128983, 20: 63025520, 21: 48129895,
         22: 51304566, 'X': 155270560,
        }
        # Randomly select a chromosome.
        chromosome = random.choice(list(chromosome_dicc.keys()))
        # Randomly generate a start and end position.
        start = np.random.randint((chromosome_dicc[chromosome] - 9_999_999))
        end = start + 10_000_000
        # Identify the variants that fall within the 10 Mb window.
        variants = np.where(((start <= Window.pos_dicc[chromosome]) & (Window.pos_dicc[chromosome] <= end)))[0]
        # If there are no variants in the 10 Mb window...
        if variants.shape[0] == 0:
            # Track the window end time.
            win_end = time.time()
            # Determine the window run time.
            win_run_time = round(win_end-win_start, 3)
            # Print the window progress report.
            sys.stdout.write(
                'window run time = {0} seconds'.format(
                    win_run_time,
                )+'\n',
            )
        # Else...
        else:
            # Determine the alternative allele frequencies for P1, P3, and P4.
            p1_freqs = Window.pop_dicc['YRI'][chromosome][variants]
            p3_freqs = Window.pop_dicc['ALT'][chromosome][variants]
            p4_freqs = Window.pop_dicc['ANC'][chromosome][variants]
            # For every P2 population....
            for p2 in p2_list:
                # Calculate site pattern counts for the bootstrapped window and store the results.
                abba, baba, baaa, abaa, abba_hom, baba_hom, baaa_hom, abaa_hom = Window.human_freq_site_patterns(
                    P1=p1_freqs,
                    P2=Window.pop_dicc[p2][chromosome][variants],
                    P3=p3_freqs,
                    P4=p4_freqs,
                )
                results_dicc[p2]['ABBA'] += abba
                results_dicc[p2]['ABBA_HOM'] += abba_hom
                results_dicc[p2]['BABA'] += baba
                results_dicc[p2]['BABA_HOM'] += baba_hom
                results_dicc[p2]['BAAA'] += baaa
                results_dicc[p2]['BAAA_HOM'] += baaa_hom
                results_dicc[p2]['ABAA'] += abaa
                results_dicc[p2]['ABAA_HOM'] += abaa_hom
            # Track the window end time.
            win_end = time.time()
            # Determine the window run time.
            win_run_time = round(win_end-win_start, 3)
            # Print the window progress report.
            sys.stdout.write(
                'window run time = {0} seconds'.format(
                    win_run_time,
                )+'\n',
            )
        return results_dicc

    # Define a function to do the work and return the window's work
    def window_worker(window):
        return window.do_work()

# Intialize the P2 populations.
tgp_p2_list = [
    'BEB', 'STU', 'ITU', 'PJL', 'GIH',
    'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
    'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
    'PEL', 'MXL', 'CLM', 'PUR',
]
# Intialize a dictionary to store bootstrapped values.
tgp_results = {}
# For all P2 populations...
for p2 in tgp_p2_list:
    # Fill the dictionary with empty site pattern lists.
    tgp_results[p2] = {
         'ABBA': 0, 'ABBA_HOM': 0,
         'BABA': 0, 'BABA_HOM': 0,
         'BAAA': 0, 'BAAA_HOM': 0,
         'ABAA': 0, 'ABAA_HOM': 0,
    }
# Build the instance of pop_dicc and pos_dicc.
Window.initialize_big_dictionaries()
# Initialize an empty list to store instances of window.
windows = []
# For 303 iterations...
for _ in range(303):
    # Create a window object and append the object list.
    windows.append(Window())
# Pool as many process as cpus.
pool = multiprocessing.Pool(multiprocessing.cpu_count())
# Assign instances of windows to pooled processes.
results = pool.map(Window.window_worker, windows)
# Clean up from multiprocessing.
pool.close()
pool.join()
# For all the result dictionaries...
for dicc in results:
    # For all P2 populations...
    for p2 in tgp_p2_list:
        tgp_results[p2]['ABBA'] += dicc[p2]['ABBA']
        tgp_results[p2]['ABBA_HOM'] += dicc[p2]['ABBA_HOM']
        tgp_results[p2]['BABA'] += dicc[p2]['BABA']
        tgp_results[p2]['BABA_HOM'] += dicc[p2]['BABA_HOM']
        tgp_results[p2]['BAAA'] += dicc[p2]['BAAA']
        tgp_results[p2]['BAAA_HOM'] += dicc[p2]['BAAA_HOM']
        tgp_results[p2]['ABAA'] += dicc[p2]['ABAA']
        tgp_results[p2]['ABAA_HOM'] += dicc[p2]['ABAA_HOM']
# Intialize a results file.
results_file = open('./bootstraps/tgp_bootstrapped_site_patterns_rep_{0}.csv'.format(str(sys.argv[1])), 'w')
# For all P2 populations...
for p2 in tgp_p2_list:
    # Compile the results for that population.
    result_line = [
        p2,
        str(tgp_results[p2]['ABBA']),
        str(tgp_results[p2]['BABA']),
        str(tgp_results[p2]['BAAA']),
        str(tgp_results[p2]['ABAA']),
        str(tgp_results[p2]['ABBA_HOM']),
        str(tgp_results[p2]['BABA_HOM']),
        str(tgp_results[p2]['BAAA_HOM']),
        str(tgp_results[p2]['ABAA_HOM']),
    ]
    # Write the results line to the results file.
    results_file.write(','.join(result_line)+'\n')
# Close the results file.
results_file.close()
