### Dependencies ###
import allel
import numpy as np
import numcodecs
import pandas as pd
import sys
import zarr
### sys.argv[1] = target population ###
### sys.argv[2] = chromosome ###


# Define subsetting function.
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
    OUTPUT: Haplodify samples for the specific pop.
    ###########################################################################
    """
    # Identify individuals in the target population.
    vals = meta_info['pop'].isin([pop]).values
    # Subset the meta info for only the target population.
    subset = meta_info[vals]
    # Get the index values of the indivuals in the target population.
    idx = subset.index.tolist()
    # Subset the genotype array for the target population and randomly sample
    # one allele.
    subset_geno_mat = geno_mat.take(idx, axis=1).haploidify_samples()
    return subset_geno_mat


# Define a function to calculate site pattern counts.
def site_patterns(
    p1,
    p2,
    p3,
):
    """
    ###########################################################################
    INPUT
        p1: P1 lineage.
        p2: P2 lineage.
        p3: P3 lineage.
    ---------------------------------------------------------------------------
    OUTPUT: Genome counts of ABBA, BABA, BBAA, BAAA, ABAA, and AABA sites.
    ###########################################################################
    """
    # Intialize site pattern counts.
    abba = 0
    baba = 0
    bbaa = 0
    baaa = 0
    abaa = 0
    aaba = 0
    # Loop through every site and calculate site pattern counts.
    for site in range(p2.shape[0]):
        if (p1[site] == 0) and (p2[site] == 1) and (p3[site] == 1):
            abba += 1
        elif (p1[site] == 1) and (p2[site] == 0) and (p3[site] == 1):
            baba += 1
        elif (p1[site] == 1) and (p2[site] == 1) and (p3[site] == 0):
            bbaa += 1
        elif (p1[site] == 1) and (p2[site] == 0) and (p3[site] == 0):
            baaa += 1
        elif (p1[site] == 0) and (p2[site] == 1) and (p3[site] == 0):
            abaa += 1
        elif (p1[site] == 0) and (p2[site] == 0) and (p3[site] == 1):
            aaba += 1
        else:
            continue
    return abba, baba, bbaa, baaa, abaa, aaba

# Define a function to calculate site pattern counts.
def calculate_human_trio_results(
    p1,
    p2,
    p3,
    p4,
):
    """
    ###########################################################################
    INPUT
        p1: P1 whole-genome sequence.
        p2: P2 whole-genome sequence.
        p3: P3 whole-genome sequence.
        p4: Outgroup whole-genome sequence.
    ---------------------------------------------------------------------------
    OUTPUT: Genome wide site pattern counts.
    ###########################################################################
    """
    # Initialize the results array for all trios.
    trio_results = np.empty((0, 8))
    # Identify the indicies where the ancestral allele is the alternative allele.
    polarize_sites_idx = np.where(p4 == 1)[0]
    # Define a lambda function to polarize sites.
    polarize = lambda x: abs(x - 1)
    # Polarize the genotypes for P1 and P3 such that the P4 is always 0.
    polarized_p1 = np.array([polarize(p1[site_idx]) if site_idx in polarize_sites_idx else p1[site_idx] for site_idx in range(len(p1))])
    polarized_p3 = np.array([polarize(p3[site_idx]) if site_idx in polarize_sites_idx else p3[site_idx] for site_idx in range(len(p3))])
    # Loop through every individual in the target population.
    for ind in range(p2.shape[1]):
        # Extract the target individual.
        p2_ind = p2[:, ind]
        # Polarize the genotypes for the target individual such that the P4
        # is always 0.
        polarized_p2 = np.array([polarize(p2_ind[site_idx]) if site_idx in polarize_sites_idx else p2_ind[site_idx] for site_idx in range(len(p2_ind))])
        # Calculate site pattern counts.
        abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
            p1=polarized_p1,
            p2=polarized_p2,
            p3=polarized_p3,
        )
        # Append the results array.
        trio_results = np.append(
            trio_results,
            [np.asarray([str(sys.argv[1]), int(sys.argv[2]), abba, baba, bbaa, baaa, abaa, aaba])], 
            axis=0,
        )
    return trio_results


# Define the zarr array path.
zarr_path = './zarr_arrays/tgp_chimp_altai_merged_filtered_biallelic_chr{0}.zarr'.format(str(sys.argv[2]))

# Read in zarr array.
callset = zarr.open_group(zarr_path, mode='r')

# Extract the genotype information.
gt_zarr = callset[str(sys.argv[2])+'/calldata/GT']

# Convert the gentoype information to a genotype array.
gt = allel.GenotypeArray(gt_zarr)

# Define the metadata path.
meta_path = './zarr_arrays/tgp_chimp_altai_info.txt'

# Read in the metadata as a pandas dataframe.
meta_data = pd.read_csv(meta_path, sep='\t', usecols=['sample', 'pop', 'super_pop'])

# Extract the genotype information for a single YRI individual.
NA18486_VAL = meta_data['sample'].isin(['NA18486']).values
NA18486_SUB = meta_data[NA18486_VAL]
NA18486_IDX = NA18486_SUB.index.tolist()
yri = gt.take(NA18486_IDX, axis=1).haploidify_samples()

# Extract the genotype information for the other focal populations.
altai = pop_subset(meta_data, gt, 'ALT')
chimp = pop_subset(meta_data, gt, 'ANC')
p2_pop = pop_subset(meta_data, gt, str(sys.argv[1]))

# Calculate the site pattern counts for every P2 individual.
site_pattern_counts = calculate_human_trio_results(yri, p2_pop, altai, chimp)

# Export the results as a csv file.
np.savetxt(
    './tgp_trios/{0}_trios_chr_{1}_site_pattern_counts.csv'.format(str(sys.argv[1]).lower(), str(sys.argv[2])),
    site_pattern_counts,
    fmt='%s',
    delimiter=',',
    newline='\n',
)
