# Import packages.
import gzip
import random
import sys
### sys.argv[1] = chromosome ###

# Define a function to build population and individual dictionaries.
def tgp_diccs(
    pop_list,
    tgp_meta_file,
):
    """
    ###########################################################################
    INPUT
        pop_list: List of populations to do site pattern counts for.
        tgp_meta_file: TGP panel file.
    ---------------------------------------------------------------------------
    OUTPUT
        pop_dicc: Dictionary of sample indicies per population.
        ind_dicc: Dictionary of site patterns per sample.
    ###########################################################################
    """
    # Intialize empty dictionaries.
    pop_dicc = {}
    ind_dicc = {}
    # For every target population.
    for pop in pop_list:
        # Fill the dictionary with the pops as keys.
        pop_dicc[pop] = []
    # Open the meta info file.
    with open(tgp_meta_file, 'r') as pop_data:
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
                # Intialize a site pattern dictionary for that indiviual.
                ind_dicc[col_counter] = {
                    'ABBA': 0,
                    'BABA': 0,
                    'BBAA': 0,
                    'BAAA': 0,
                    'ABAA': 0,
                    'AABA': 0,
                }
                # Move the colum counter forward.
                col_counter += 1
            # Else...
            else:
                # Move the colum counter forward.
                col_counter += 1
    return pop_dicc, ind_dicc

# Define a site pattern function.
def tgp_trios_site_patterns(
    pop_list,
    pop_dicc,
    ind_dicc,
    vcf,
):
    """
    ###########################################################################
    INPUT
        pop_list: List of populations to do site pattern counts for.
        pop_dicc: Dictionary of sample indicies per population.
        ind_dicc: Dictionary of site patterns per sample.
        vcf: tgp_chimp_{archaic}_merged_filtered_biallelic_chr{#}.vcf.gz
    ---------------------------------------------------------------------------
    OUTPUT
        pop_dicc: Dictionary of sample indicies per population with site pattern
                  counts.
        header_info: Header line split by tabs.
    ###########################################################################
    """    
    # Using the gzip package open the vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Intialize RNG.
        rng_vals = [0, 2]
        # Loop through every line in the vcf file.
        for line in data:
            # If the current line is a part of the meta info...
            if line.startswith('##'):
                # Continue to the next line in the vcf file.
                continue
            # Else-if the current line is the header line.
            elif line.startswith('#'):
                # Split the line.
                spline = line.split()
                # Save the header line.
                header_info = spline
            # Else...
            else:
                # Split the line.
                spline = line.split()
                # Randomly select a chromosome to sample.
                random_chr = random.choice(rng_vals)
                # Randomly sample a chromosome for the P1, P3, and P4 indiviuals.
                p1 = spline[1764][random_chr] # YRI individual.
                p3 = spline[-2][random_chr] # Archaic individual.
                p4 = spline[-1][random_chr] # EPO ancestral allele call.
                # For every target population...
                for pop in pop_list:
                    # For every individual in the target population.
                    for ind in pop_dicc[pop]:
                        # Randomly sample a chromosome.
                        p2 = spline[ind][random_chr]
                        # Determine the site pattern for that individual.
                        # ABBA.
                        if ((p1 == '0') and (p2 == '1') and (p3 == '1') and (p4 == '0')):
                            ind_dicc[ind]['ABBA'] += 1
                        elif ((p1 == '1') and (p2 == '0') and (p3 == '0') and (p4 == '1')):
                            ind_dicc[ind]['ABBA'] += 1
                        # BABA.
                        elif ((p1 == '1') and (p2 == '0') and (p3 == '1') and (p4 == '0')):
                            ind_dicc[ind]['BABA'] += 1
                        elif ((p1 == '0') and (p2 == '1') and (p3 == '0') and (p4 == '1')):
                            ind_dicc[ind]['BABA'] += 1
                        # BBAA.
                        elif ((p1 == '1') and (p2 == '1') and (p3 == '0') and (p4 == '0')):
                            ind_dicc[ind]['BBAA'] += 1
                        elif ((p1 == '0') and (p2 == '0') and (p3 == '1') and (p4 == '1')):
                            ind_dicc[ind]['BBAA'] += 1
                        # BAAA.
                        elif ((p1 == '1') and (p2 == '0') and (p3 == '0') and (p4 == '0')):
                            ind_dicc[ind]['BAAA'] += 1
                        elif ((p1 == '0') and (p2 == '1') and (p3 == '1') and (p4 == '1')):
                            ind_dicc[ind]['BAAA'] += 1
                        # ABAA.
                        elif ((p1 == '0') and (p2 == '1') and (p3 == '0') and (p4 == '0')):
                            ind_dicc[ind]['ABAA'] += 1
                        elif ((p1 == '1') and (p2 == '0') and (p3 == '1') and (p4 == '1')):
                            ind_dicc[ind]['ABAA'] += 1
                        # AABA.
                        elif ((p1 == '0') and (p2 == '0') and (p3 == '1') and (p4 == '0')):
                            ind_dicc[ind]['AABA'] += 1
                        elif ((p1 == '1') and (p2 == '1') and (p3 == '0') and (p4 == '1')):
                            ind_dicc[ind]['AABA'] += 1
                        else:
                            continue
    return ind_dicc, header_info

# Load the sys args.
chrom = str(sys.argv[1])

# Load in the vcf file.
chrom_vcf = './tgp_altai_merged_filtered_biallelic_chr{0}.vcf.gz'.format(chrom)

# Load in the meta info.
tgp_meta_info = './tgp_altai_ancestor_info.txt'

# Define target populations.
target_pops = [
    'CEU', 'FIN', 'GBR', 'IBS', 'TSI',
    'CHB', 'CHS', 'CDX', 'JPT', 'KHV',
    'BEB', 'GIH', 'ITU', 'PJL', 'STU',
    'CLM', 'MXL', 'PEL', 'PUR',
]

# Create dictionaries.
tgp_pop_dicc, tgp_ind_dicc = tgp_diccs(target_pops, tgp_meta_info)

# Calculate site patterns per trio.
tgp_trios_dicc, tgp_header = tgp_trios_site_patterns(
    target_pops,
    tgp_pop_dicc,
    tgp_ind_dicc,
    chrom_vcf,
)

# Intialize a results file.
results_file = open('./tgp_trios/tgp_altai_trios_chr_{0}_site_pattern_counts.csv'.format(chrom), 'w')
# For all target populations...
for pop in target_pops:
    # For all individuals in that population...
    for ind in tgp_pop_dicc[pop]:
        # Compile the results for that individual.
        result_line = [
            pop,
            tgp_header[ind],
            chrom,
            str(tgp_trios_dicc[ind]['ABBA']),
            str(tgp_trios_dicc[ind]['BABA']),
            str(tgp_trios_dicc[ind]['BBAA']),
            str(tgp_trios_dicc[ind]['BAAA']),
            str(tgp_trios_dicc[ind]['ABAA']),
            str(tgp_trios_dicc[ind]['AABA']),
        ]
        # Write the results line to the results file.
        results_file.write(','.join(result_line)+'\n')
# Close the results file.
results_file.close()
