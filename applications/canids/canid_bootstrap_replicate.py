# Import packages.
import gzip
import numpy as np
import sys


# Define a function to preform one bootstrap replicate.
def canid_bootstrap_replicate(
    p1, 
    p2,
    p3,
):
    """
    ###########################################################################
    INPUT:
        p1: ID of the P1 individual.
        p2: ID of the P2 individual.
        p3: ID of the P3 individual.
    ---------------------------------------------------------------------------
    OUTPUT: List of site pattern counts from one bootstrapped replicate.
    ###########################################################################
    """
    # Intialize a dictionary with sample indicies.
    sample_dicc = {
        'BSJ': 9, 'DNG': 10, 'ISW': 11,
        'CRW': 12, 'CHW': 13, 'GLJ': 14,
    }
    # Determine sample indicies.
    p1_idx = sample_dicc[p1]
    p2_idx = sample_dicc[p2]
    p3_idx = sample_dicc[p3]
    p4_idx = sample_dicc['GLJ']
    # Intialize a site pattern dictionarry.
    site_pattern_dicc = {
        'ABBA': 0, 'BABA': 0, 'BBAA': 0,
        'BAAA': 0, 'ABAA': 0, 'AABA': 0,
    }
    # Intialize a dictionary of chromosome lengths for CanFam3.1 assembly.
    chromosome_dicc = {
        1: 122678785, 2: 85426708, 3: 91889043, 4: 88276631, 5: 88915250,
        6: 77573801, 7: 80974532, 8: 74330416, 9: 61074082, 10: 69331447,
        11: 74389097, 12: 72498081, 13: 63241923, 14: 60966679, 15: 64190966,
        16: 59632846, 17: 64289059, 18: 55844845, 19: 53741614, 20: 58134056,
        21: 50858623, 22: 61439934, 23: 52294480, 24: 47698779, 25: 51628933,
        26: 38964690, 27: 45876710, 28: 41182112, 29: 41845238, 30: 40214260,
        31: 39895921, 32: 38810281, 33: 31377067, 34: 42124431, 35: 26524999,
        36: 30810995, 37: 30902991, 38: 23914537,
    }
    # Intialize empty dictionary for variable positions.
    pos_dicc = {}
    # Fill the positions dictionary with positions arrays.
    for chrom in list(chromosome_dicc.keys()):
        pos_dicc[chrom] = np.loadtxt(
            './canid_positions_chr{0}.csv.gz'.format(chrom),
            dtype=int, delimiter=',',
        )
    # For each bootstrap window (220 x 10 Mb)...
    for j in range(220):
        # Randomly select a chromosome.
        chromosome = np.random.choice(list(chromosome_dicc.keys()))
        # Randomly generate a start and end position.
        start = np.random.randint((chromosome_dicc[chromosome] - 9_999_999))
        end = start + 10_000_000
        # Identify the indicies for this window.
        variants = np.where(((start <= pos_dicc[chromosome]) & (pos_dicc[chromosome] <= end)))[0]
        # If there are no variants in the window...
        if variants.shape[0] == 0:
            # Continue to the next window.
            continue
        # Else...
        else:
            # Using the gzip package open the vcf file.
            vcf = './canid_merged_filtered_chr{0}.vcf.gz'.format(chromosome)
            with gzip.open(vcf, 'rt') as data:
                # Intialize rng values.
                rng_vals = [0, 2]
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
                        # Grab the sample info.
                        p1_info = spline[p1_idx]
                        p2_info = spline[p2_idx]
                        p3_info = spline[p3_idx]
                        p4_info = spline[p4_idx]
                        # If any sample has missing data...
                        if '.:.:.' in [p1_info, p2_info, p3_info, p4_info]:
                            # Move the sites index counter forward.
                            idx_counter += 1
                            # Continue to the next line.
                            continue
                        # Else...
                        else:
                            # Grab the sample filter info.
                            p1_sf = int(p1_info[-1])
                            p2_sf = int(p2_info[-1])
                            p3_sf = int(p3_info[-1])
                            p4_sf = int(p4_info[-1])
                            # If all samples pass the filtering...
                            if (p1_sf+p2_sf+p3_sf+p4_sf) == 4:
                                # Move the sites index counter forward.
                                idx_counter += 1
                                # Randomly select one chromosome to sample.
                                random_chrom = np.random.choice(rng_vals)
                                # Randomly sample one chromosome for each sample.
                                p1_chrom = p1_info[random_chrom]
                                p2_chrom = p2_info[random_chrom]
                                p3_chrom = p3_info[random_chrom]
                                p4_chrom = p4_info[random_chrom]
                                # Determine the site pattern.
                                # ABBA.
                                if ((p1_chrom == '0') and (p2_chrom == '1') and (p3_chrom == '1') and (p4_chrom == '0')):
                                    site_pattern_dicc['ABBA'] += 1
                                elif ((p1_chrom == '1') and (p2_chrom == '0') and (p3_chrom == '0') and (p4_chrom == '1')):
                                    site_pattern_dicc['ABBA'] += 1
                                # BABA.
                                elif ((p1_chrom == '1') and (p2_chrom == '0') and (p3_chrom == '1') and (p4_chrom == '0')):
                                    site_pattern_dicc['BABA'] += 1
                                elif ((p1_chrom == '0') and (p2_chrom == '1') and (p3_chrom == '0') and (p4_chrom == '1')):
                                    site_pattern_dicc['BABA'] += 1
                                # BBAA.
                                elif ((p1_chrom == '1') and (p2_chrom == '1') and (p3_chrom == '0') and (p4_chrom == '0')):
                                    site_pattern_dicc['BBAA'] += 1
                                elif ((p1_chrom == '0') and (p2_chrom == '0') and (p3_chrom == '1') and (p4_chrom == '1')):
                                    site_pattern_dicc['BBAA'] += 1
                                # BAAA.
                                elif ((p1_chrom == '1') and (p2_chrom == '0') and (p3_chrom == '0') and (p4_chrom == '0')):
                                    site_pattern_dicc['BAAA'] += 1
                                elif ((p1_chrom == '0') and (p2_chrom == '1') and (p3_chrom == '1') and (p4_chrom == '1')):
                                    site_pattern_dicc['BAAA'] += 1
                                # ABAA.
                                elif ((p1_chrom == '0') and (p2_chrom == '1') and (p3_chrom == '0') and (p4_chrom == '0')):
                                    site_pattern_dicc['ABAA'] += 1
                                elif ((p1_chrom == '1') and (p2_chrom == '0') and (p3_chrom == '1') and (p4_chrom == '1')):
                                    site_pattern_dicc['ABAA'] += 1
                                # AABA.
                                elif ((p1_chrom == '0') and (p2_chrom == '0') and (p3_chrom == '1') and (p4_chrom == '0')):
                                    site_pattern_dicc['AABA'] += 1
                                elif ((p1_chrom == '1') and (p2_chrom == '1') and (p3_chrom == '0') and (p4_chrom == '1')):
                                    site_pattern_dicc['AABA'] += 1
                                # Else...
                                else:
                                    # Continue to the next line.
                                    continue
                            # Else...
                            else:
                                # Move the sites index counter forward.
                                idx_counter += 1
                                # Continue to the next line.
                                continue
    return site_pattern_dicc


# Perform one bootstrap replicate.
dng_bsj_isw_bs_dicc = canid_bootstrap_replicate(p1='DNG', p2='BSJ', p3='ISW')

# Intialize a results file.
results_file = open('./bootstraps/canid_bootstrapped_site_patterns_rep_{0}.csv'.format(str(sys.argv[1])), 'w')

# Compile the results as a list.
dng_bsj_isw_line = [
    str(dng_bsj_isw_bs_dicc['ABBA']), str(dng_bsj_isw_bs_dicc['BABA']), str(dng_bsj_isw_bs_dicc['BBAA']),
    str(dng_bsj_isw_bs_dicc['BAAA']), str(dng_bsj_isw_bs_dicc['ABAA']), str(dng_bsj_isw_bs_dicc['AABA']),
]

# Write the results line to the results file.
results_file.write(','.join(dng_bsj_isw_line)+'\n')

# Close the results file.
results_file.close()