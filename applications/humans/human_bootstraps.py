### Dependencies ###
import gzip
import numpy as np
import random
import sys
import time
### sys.argv[1] = bootstrap replicate ID ###

# Define a function to preform one bootstrap replicate.
def human_bootstrap_replicate():
    """
    ###########################################################################
    INPUT: N/A
    ---------------------------------------------------------------------------
    OUTPUT: List of site pattern counts from one bootstrapped replicate.
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
    # Intialize the P2 populations.
    p2_list = [
        'BEB', 'STU', 'ITU', 'PJL', 'GIH',
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
        'PEL', 'MXL', 'CLM', 'PUR',
    ]
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
    # Intialize a dictionary for sample indicies.
    pop_dicc = {}
    # For every target population.
    for pop in pop_list:
        # Fill the dictionary with the pops as keys.
        pop_dicc[pop] = []
    # Open the meta info file.
    with open('./zarr_arrays/tgp_chimp_altai_info.txt', 'r') as pop_data:
        # Skip the header line.
        next(pop_data)
        # Intialize a column counter.
        col_counter = 9
        # For every line...
        for line in pop_data:
            # Split the line.
            spline = line.split()
            # Grab the current population.
            pop = spline[1]
            # If the population is a target population...
            if pop in pop_list:
                # Append the individual to the population list.
                pop_dicc[pop].append(col_counter)
                # Move the colum counter forward.
                col_counter += 1
            # Else...
            else:
                # Move the colum counter forward.
                col_counter += 1
    # Intialize a dictionary to store bootstrapped values.
    results_dicc = {}
    # For all P2 populations.
    for p2 in p2_list:
        # Fill the dictionary by intializing site patterns.
        results_dicc[p2] = {
            'ABBA': 0, 'ABBA_HOM': 0,
            'BABA': 0, 'BABA_HOM': 0,
            'BAAA': 0, 'BAAA_HOM': 0,
            'ABAA': 0, 'ABAA_HOM': 0,
        }
    # Intialize empty dictionary for variable positions.
    pos_dicc = {}
    # Fill the positions dictionary with positions arrays.
    for chrom in list(chromosome_dicc.keys()):
        pos_dicc[chrom] = np.loadtxt(
            './zarr_arrays/tgp_chimp_altai_merged_filtered_biallelic_chr{0}_positions.csv.gz'.format(chrom),
            dtype=int, delimiter=',',
        )
    # Track the total run time.
    total_run_time = 0
    # For each bootstrap window (303 x 10 Mb)...
    for j in range(303):
        # Track the window start time.
        win_start = time.time()
        # Randomly select a chromosome.
        chromosome = random.choice(list(chromosome_dicc.keys()))
        # Randomly generate a start and end position.
        start = np.random.randint((chromosome_dicc[chromosome] - 9_999_999))
        end = start + 10_000_000
        # Print the window information.
        sys.stdout.write(
            'window: {0}; chr: {1}; start: {2}; end: {3}'.format(
                j, chromosome, start, end,
            )+'\n',
        )
        # Identify the indicies for this window.
        variants = np.where(((start <= pos_dicc[chromosome]) & (pos_dicc[chromosome] <= end)))[0]
        # If there are no variants in the window...
        if variants.shape[0] == 0:
            # Track the window end time.
            win_end = time.time()
            # Determine the window run time.
            win_run_time = round(win_end-win_start, 3)
            # Print the window progress report.
            sys.stdout.write(
                'window = {0}; run time = {1} seconds'.format(
                    j, round(win_run_time, 3),
                )+'\n',
            )
            # Keep track of the total run time
            total_run_time += win_run_time
            # Continue to the next window.
            continue
        # Else...
        else:
            # Using the gzip package open the vcf file.
            vcf = './zarr_arrays/tgp_chimp_altai_merged_filtered_biallelic_chr{0}.vcf.gz'.format(chromosome)
            with gzip.open(vcf, 'rt') as data:
                # Intialize a sites index counter.
                idx_counter = 0
                # For every line in the vcf file...
                for line in data:
                    # If the current line is a part of the meta info or header...
                    if line.startswith('#'):
                        # Continue to the next line in the vcf file.
                        continue
                    # Else-if the current site is larger than the end of the window...
                    elif idx_counter > variants[-1]:
                        # Break and continue to the next window.
                        break
                    # Else-if the current site is smaller than the beginning of the window...
                    elif idx_counter < variants[0]:
                        # Move the sites index counter forward.
                        idx_counter += 1
                        # Continue to the next line in the vcf file.
                        continue
                    # Else...
                    else:
                        # Split the line.
                        spline = line.split()
                        # Intialize alternative allele counters for P1.
                        p1_alt_count = 0
                        # Intialize the number of sampled chromosomes for P1.
                        p1_chroms = len(pop_dicc['YRI']) * 2
                        # For all individuals in P1...
                        for yri in pop_dicc['YRI']:
                            # Update the count of alternative alleles.
                            p1_alt_count += spline[yri][0:3].count('1')
                        # Determine the alternative allele frequency for P1.
                        p1_freq = p1_alt_count / p1_chroms   
                        # Determine the alternative allele frequency for P3.   
                        p3_freq = (spline[pop_dicc['ALT'][0]][0:3].count('1')) / 2    
                        # Determine the alternative allele frequency for P4.   
                        p4_freq = (spline[pop_dicc['ANC'][0]][0:3].count('1')) / 2      
                        # For every P2 population...
                        for p2 in p2_list:
                            # Intialize the alternative allele counter.
                            p2_alt_count = 0
                            # Intialize the number of samples.
                            p2_chroms = len(pop_dicc[p2]) * 2
                            # For every individual...
                            for ind in pop_dicc[p2]:
                                # Update the count of alternative alleles.
                                p2_alt_count += spline[ind][0:3].count('1')
                            # Determine the alternative allele frequency.
                            p2_freq = p2_alt_count / p2_chroms
                            # If the ancestral allele is the refernce allele.
                            if p4_freq == 0.0:
                                # Calculate site patterns.
                                results_dicc[p2]['ABBA'] += (1 - p1_freq) * p2_freq * p3_freq
                                results_dicc[p2]['ABBA_HOM'] += (1 - p1_freq) * p3_freq * p3_freq
                                results_dicc[p2]['BABA'] += p1_freq * (1 - p2_freq) * p3_freq
                                results_dicc[p2]['BABA_HOM'] += p1_freq * (1 - p3_freq) * p3_freq
                                results_dicc[p2]['BAAA'] += p1_freq * (1 - p2_freq) * (1 - p3_freq)
                                results_dicc[p2]['BAAA_HOM'] += p1_freq * (1 - p3_freq) * (1 - p3_freq)
                                results_dicc[p2]['ABAA'] += (1 - p1_freq) * p2_freq * (1 - p3_freq)
                                results_dicc[p2]['ABAA_HOM'] += (1 - p1_freq) * p3_freq * (1 - p3_freq)
                            # Else-if the ancestral allele is the alternative allele.
                            elif p4_freq == 1.0:
                                # Calculate site patterns.
                                results_dicc[p2]['ABBA'] += p1_freq * (1 - p2_freq) * (1 - p3_freq)
                                results_dicc[p2]['ABBA_HOM'] += p1_freq * (1 - p3_freq) * (1 - p3_freq)
                                results_dicc[p2]['BABA'] += (1 - p1_freq) * p2_freq * (1 - p3_freq)
                                results_dicc[p2]['BABA_HOM'] += (1 - p1_freq) * p3_freq * (1 - p3_freq)
                                results_dicc[p2]['BAAA'] += (1 - p1_freq) * p2_freq * p3_freq
                                results_dicc[p2]['BAAA_HOM'] += (1 - p1_freq) * p3_freq * p3_freq
                                results_dicc[p2]['ABAA'] += p1_freq * (1 - p2_freq) * p3_freq
                                results_dicc[p2]['ABAA_HOM'] += p1_freq * (1 - p3_freq) * p3_freq
                    # Move the sites index counter forward.
                    idx_counter += 1
        # Track the window end time.
        win_end = time.time()
        # Determine the window run time.
        win_run_time = round(win_end-win_start, 3)
        # Print the window progress report.
        sys.stdout.write(
            'window = {0}; run time = {1} minutes'.format(
                j, round(win_run_time/60, 3),
            )+'\n',
        )
        # Keep track of the total run time.
        total_run_time += win_run_time
    # Print the replicate report.
    sys.stdout.write(
        'completed bootstrapped replicate in {0} minutes'.format(
            round(total_run_time/60, 3),
        )+'\n',
    )
    return results_dicc


# Run one bootstrapped replicate.
bs_dicc = human_bootstrap_replicate()

# Intialize the P2 list for the results.
p2_pops = [
        'BEB', 'STU', 'ITU', 'PJL', 'GIH',
        'TSI', 'CEU', 'IBS', 'GBR', 'FIN',
        'CHB', 'KHV', 'CHS', 'JPT', 'CDX',
        'PEL', 'MXL', 'CLM', 'PUR',
    ]

# Intialize a results file.
results_file = open('./bootstraps/tgp_bootstrapped_site_patterns_rep_{0}.csv'.format(str(sys.argv[1])), 'w')
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