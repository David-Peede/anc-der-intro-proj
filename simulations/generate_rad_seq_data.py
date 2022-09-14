### Dependencies ###
from functools import reduce
import msprime
import numpy as np
import sys
### sys.argv[1] = the admixture proportion ###
### sys.argv[2] = replicate ID ###


# Define derived allele frequency function.
def derived_allele_freq(
    genotype_matrix,
    taxon,
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genotypes from msprime.
        taxon: Lineage to calculate the derived allele frequencies for.
    ---------------------------------------------------------------------------
    OUTPUT: Derived allele frequency array for the lineage of interest.
    ---------------------------------------------------------------------------
    NOTE: This function will return counts when the sample size is 1.
    ###########################################################################
    """
    # Calculate derived allele frequencies.
    freq_array = genotype_matrix[:, taxon].sum(axis=1)/float(len(taxon))
    return freq_array

# Define site pattern function.
def site_patterns(
    genotype_matrix,
    p1_idx=[0],
    p2_idx=[1],
    p3_idx=[2],
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genotypes from msprime.
        p1_idx: P1 lineage.
        p2_idx: P2 lineage.
        p3_idx: P3 lineage.
    ---------------------------------------------------------------------------
    OUTPUT: Genome counts of ABBA, BABA, BBAA, BAAA, ABAA, and AABA sites.
    ###########################################################################
    """
    # Calculate derived allele frequencies per taxa.
    p1 = derived_allele_freq(genotype_matrix, p1_idx)
    p2 = derived_allele_freq(genotype_matrix, p2_idx)
    p3 = derived_allele_freq(genotype_matrix, p3_idx)
    # Calculate site pattern counts.
    abba = ((1 - p1) * (p2) * (p3)).sum()
    baba = ((p1) * (1 - p2) * (p3)).sum()
    bbaa = ((p1) * (p2) * (1 - p3)).sum()
    baaa = ((p1) * (1 - p2) * (1 - p3)).sum()
    abaa = ((1 - p1) * (p2) * (1 - p3)).sum()
    aaba = ((1 - p1) * (1 - p2) * (p3)).sum()
    return abba, baba, bbaa, baaa, abaa, aaba

# Define a function to conduct the RADseq experiments.
def rad_loci_bootstrapping(
    genotype_matrix,
    variable_positions,
    n_replicates,
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genome.
        variable_positions: Variable positions associated with the simulated
                            genome.
        n_replicates: Number of bootstrap replicates.
    ---------------------------------------------------------------------------
    OUTPUT: RADseq results from bootstrapped replicates.
    ###########################################################################
    """
    # Intialize arrays to store bootstrapped values.
    abba_array = np.array([])
    baba_array = np.array([])
    bbaa_array = np.array([])
    baaa_array = np.array([])
    abaa_array = np.array([])
    aaba_array = np.array([])
    # Intialize lists for the cut sequence and number rad sites.
    cut_seq_list = []
    n_rad_sites_list = []
    # For each bootstrap replicate.
    for i in range(n_replicates):
        # Randomly generate a cutting sequence.
        cut_seq = ''.join(np.random.binomial(1, 0.5, 12).astype('str').tolist())
        # Intialize variable positions sequence.
        p1_var_seq = ''.join(genotype_matrix[:, 0].astype('str').tolist())
        p2_var_seq = ''.join(genotype_matrix[:, 1].astype('str').tolist())
        p3_var_seq = ''.join(genotype_matrix[:, 2].astype('str').tolist())
        # Intialize empty lists to store cut sites indicies.
        p1_cut_sites_idx = []
        p2_cut_sites_idx = []
        p3_cut_sites_idx = []
        # Intialize the cut lengths.
        cut_length = 12
        # For every variable position...
        for pos in range(genotype_matrix.shape[0] - cut_length):
            # Grab the site information per population.
            p1_sites = p1_var_seq[pos:pos+cut_length]
            p2_sites = p2_var_seq[pos:pos+cut_length]
            p3_sites = p3_var_seq[pos:pos+cut_length]
            # If the P1 individual has the cut sequence...
            if p1_sites == cut_seq:
                # Append the current position.
                p1_cut_sites_idx.append(pos)
            # If the P2 individual has the cut sequence...
            if p2_sites == cut_seq:
                # Append the current position.
                p2_cut_sites_idx.append(pos)
            # If the P3 individual has the cut sequence...
            if p3_sites == cut_seq:
                # Append the current position.
                p3_cut_sites_idx.append(pos)
        # Find all of cut site indicies among all three populations.
        all_cut_sites_idx = reduce(
            np.union1d,
            (np.array(p1_cut_sites_idx), np.array(p2_cut_sites_idx), np.array(p3_cut_sites_idx)),
        )
        # If there are cut sites...
        if all_cut_sites_idx.size > 0:
            # Find all cut site positions.
            all_cut_positions = variable_positions[all_cut_sites_idx.astype(int)]
            # Find all of the indicies where the cut sites are within 100 bp of the next cut site.
            within_100bp_idx = np.where(np.diff(all_cut_positions, append=99_999_999) <= 100)[0]
            # If there are cut sites to prune...
            if within_100bp_idx.size > 0:
                # If the last varible position is a cut site and needs to be pruned...
                if (within_100bp_idx.size - 1) in within_100bp_idx:
                    # Get the indicies of the next cut site, ie cut_positions[idx+1].
                    next_100bp_idx = within_100bp_idx[:-1] + 1
                    # Combine the original and next cut site indicies to be pruned.
                    prune_idx = np.union1d(within_100bp_idx, next_100bp_idx)
                    # Filter the cut sites indicies.
                    filtered_cut_sites_idx = np.setdiff1d(all_cut_sites_idx, prune_idx)
                # Else...
                else:
                    # Get the positions of the next cut site, ie cut_positions[idx+1].
                    next_100bp_idx = within_100bp_idx + 1
                    # Combine the original and next cut site indicies to be pruned.
                    prune_idx = np.union1d(within_100bp_idx, next_100bp_idx)
                    # Filter the cut sites indicies.
                    filtered_cut_sites_idx = np.setdiff1d(all_cut_sites_idx, prune_idx)
            # Else...
            else:
                # Make the filtered cut sites indicies the same as all the cut site indicies.
                filtered_cut_sites_idx = all_cut_sites_idx.astype(int)
            # Filter the genetoype matrix by cut sites.
            rad_genotype_matrix = genotype_matrix[filtered_cut_sites_idx.astype(int), :]
            # Calculate site pattern counts for the bootstrapped replicate and store the results.
            abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
                rad_genotype_matrix,
                p1_idx=[0],
                p2_idx=[1],
                p3_idx=[2],
            )
            abba_array = np.append(abba_array, abba)
            baba_array = np.append(baba_array, baba)
            bbaa_array = np.append(bbaa_array, bbaa)
            baaa_array = np.append(baaa_array, baaa)
            abaa_array = np.append(abaa_array, abaa)
            aaba_array = np.append(aaba_array, aaba)
            # Append the cut sequence and number of rad sites results.
            cut_seq_list.append(cut_seq)
            n_rad_sites_list.append(rad_genotype_matrix.shape[0])
        # Else...
        else:
            # Append the site pattern results with NaNs.
            abba_array = np.append(abba_array, np.nan)
            baba_array = np.append(baba_array, np.nan)
            bbaa_array = np.append(bbaa_array, np.nan)
            baaa_array = np.append(baaa_array, np.nan)
            abaa_array = np.append(abaa_array, np.nan)
            aaba_array = np.append(aaba_array, np.nan)
            # Append the cut sequence and number of rad sites results.
            cut_seq_list.append(cut_seq)
            n_rad_sites_list.append(0)
    return abba_array, baba_array, bbaa_array, baaa_array, abaa_array, aaba_array, cut_seq_list, np.array(n_rad_sites_list)


# Parse comand-line arguments.
f         = float(sys.argv[1])
rep_id    = int(sys.argv[2])

# Load the genotype matrix and variable positions.
simulated_genotype_matrix = np.loadtxt(
    './sim_outputs/{0}/geno_mats/rep_id_{1}_geno_mat.csv.gz'.format(f, rep_id),
    dtype=int, delimiter=',',
)
simulated_variable_positions = np.loadtxt(
    './sim_outputs/{0}/var_pos/rep_id_{1}_var_pos.csv.gz'.format(f, rep_id),
    delimiter=',',
)

# Perform RADseq experiments.
rad_abba, rad_baba, rad_bbaa, rad_baaa, rad_abaa, rad_aaba, rad_cut_seqs, rad_seg_sites = rad_loci_bootstrapping(
    genotype_matrix=simulated_genotype_matrix,
    variable_positions=simulated_variable_positions,
    n_replicates=1_000,
)

# Save each RADseq bootstrapped distribution as a gzipped csv.
results_prefix = './sim_outputs/{0}/rad_seq/rep_id_{1}_'.format(f, rep_id)

# Save the cut site sequences as a csv.
cut_seqs_file = open(results_prefix+'cut_seqs.csv', 'w')
cut_seqs_file.write(','.join(rad_cut_seqs))
cut_seqs_file.close()

np.savetxt(
    results_prefix+'abba.csv.gz',
    [rad_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baba.csv.gz',
    [rad_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'bbaa.csv.gz',
    [rad_bbaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'baaa.csv.gz',
    [rad_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'abaa.csv.gz',
    [rad_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'aaba.csv.gz',
    [rad_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'seg_sites.csv.gz',
    [rad_seg_sites],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)