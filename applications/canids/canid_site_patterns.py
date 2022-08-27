# Import packages.
import gzip
import numpy as np


# Define a function to calculate site patterns.
def canid_site_patterns(
    p1,
    p2,
    p3
):
    """
    ###########################################################################
    INPUT:
        p1: ID of the P1 individual.
        p2: ID of the P2 individual.
        p3: ID of the P3 individual.
    ---------------------------------------------------------------------------
    OUTPUT
        site_pattern_dicc: Dictionary of genome wide site pattern counts.
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
    # For every chromosome...
    for chrom in range(1, 39):
        # Grab the vcf file.
        vcf = './canid_merged_filtered_chr{0}.vcf.gz'.format(chrom)
        # Using the gzip package open the vcf file.
        with gzip.open(vcf, 'rt') as data:
            # Intialize rng values.
            rng_vals = [0, 2]
            # Loop through every line in the vcf file.
            for line in data:
                # If the current line is a part of the meta info...
                if line.startswith('#'):
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
                            # Continue to the next line.
                            continue
    return site_pattern_dicc


# Calculate site patterns.
dng_bsj_isw_dicc = canid_site_patterns(p1='DNG', p2='BSJ', p3='ISW')

# Intialize a results file.
results_file = open('./canid_site_pattern_counts.csv', 'w')

# Compile the results as a list.
dng_bsj_isw_line = [
    str(dng_bsj_isw_dicc['ABBA']), str(dng_bsj_isw_dicc['BABA']), str(dng_bsj_isw_dicc['BBAA']),
    str(dng_bsj_isw_dicc['BAAA']), str(dng_bsj_isw_dicc['ABAA']), str(dng_bsj_isw_dicc['AABA']),
]

# Write the results line to the results file.
results_file.write(','.join(dng_bsj_isw_line))

# Close the results file.
results_file.close()