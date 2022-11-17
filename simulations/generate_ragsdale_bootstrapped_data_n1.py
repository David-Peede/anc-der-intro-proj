### Dependencies ###
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

# Define bootstrapping function.
def ooa_arc_intro_genotype_matrix_bootstrap(
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
    OUTPUT: List of site pattern counts from bootstrapped replicates.
    ###########################################################################
    """
    # Intialize arrays to store bootstrapped values.
    # CEU.
    ceu_abba_array = np.array([])
    ceu_baba_array = np.array([])
    ceu_bbaa_array = np.array([])
    ceu_baaa_array = np.array([])
    ceu_abaa_array = np.array([])
    ceu_aaba_array = np.array([])
    # CHB.
    chb_abba_array = np.array([])
    chb_baba_array = np.array([])
    chb_bbaa_array = np.array([])
    chb_baaa_array = np.array([])
    chb_abaa_array = np.array([])
    chb_aaba_array = np.array([])
    # For each bootstrap replicate.
    for i in range(n_replicates):
        # Intialize an empty array to store the bootstrapped genome.
        bootstrapped_genome = np.empty((0, 4))
        # For each bootstrap window (1000 x 100 kb).
        for j in range(1000):
            # Randomly generate a start and end position.
            start = np.random.randint(99_900_001)
            end = start + 100_000
            # Identify the variants that fall within the 100 kb window.
            variants = np.where(((start <= variable_positions) & (variable_positions <= end)))[0]
            # Subset the 100 kb window from the genome-wide genotype matrix.
            windowed_matrix = genotype_matrix[variants, :]
            # Concatenate 1000 windows of 100 kb to make a bootstrapped genome.
            bootstrapped_genome = np.concatenate((bootstrapped_genome, windowed_matrix))
        # Calculate site pattern counts for the bootstrapped replicate and store the results.
        # CEU.
        ceu_abba, ceu_baba, ceu_bbaa, ceu_baaa, ceu_abaa, ceu_aaba = site_patterns(
            genotype_matrix=bootstrapped_genome,
            p1_idx=[0],
            p2_idx=[1],
            p3_idx=[3],
        )
        ceu_abba_array = np.append(ceu_abba_array, ceu_abba)
        ceu_baba_array = np.append(ceu_baba_array, ceu_baba)
        ceu_bbaa_array = np.append(ceu_bbaa_array, ceu_bbaa)
        ceu_baaa_array = np.append(ceu_baaa_array, ceu_baaa)
        ceu_abaa_array = np.append(ceu_abaa_array, ceu_abaa)
        ceu_aaba_array = np.append(ceu_aaba_array, ceu_aaba)
        # CHB.
        chb_abba, chb_baba, chb_bbaa, chb_baaa, chb_abaa, chb_aaba = site_patterns(
            genotype_matrix=bootstrapped_genome,
            p1_idx=[0],
            p2_idx=[2],
            p3_idx=[3],
        )
        chb_abba_array = np.append(chb_abba_array, chb_abba)
        chb_baba_array = np.append(chb_baba_array, chb_baba)
        chb_bbaa_array = np.append(chb_bbaa_array, chb_bbaa)
        chb_baaa_array = np.append(chb_baaa_array, chb_baaa)
        chb_abaa_array = np.append(chb_abaa_array, chb_abaa)
        chb_aaba_array = np.append(chb_aaba_array, chb_aaba)
    # Aggregate results.
    ceu_results = [ceu_abba_array, ceu_baba_array, ceu_bbaa_array, ceu_baaa_array, ceu_abaa_array, ceu_aaba_array]
    chb_results = [chb_abba_array, chb_baba_array, chb_bbaa_array, chb_baaa_array, chb_abaa_array, chb_aaba_array]
    return ceu_results, chb_results


# Parse comand-line arguments.
f         = float(sys.argv[1])
rep_id    = int(sys.argv[2])

# Load in the replicate data.
simulated_genotype_matrix = np.loadtxt(
    './non_iua_models/ragsdale_2019/n_1/{0}/geno_mats/rep_id_{1}_geno_mat.csv.gz'.format(f, rep_id),
    dtype=int, delimiter=',',
)
simulated_variable_positions = np.loadtxt(
    './non_iua_models/ragsdale_2019/n_1/{0}/var_pos/rep_id_{1}_var_pos.csv.gz'.format(f, rep_id),
    delimiter=',',
)

# Perform bootstrapping.
bs_ceu, bs_chb, = ooa_arc_intro_genotype_matrix_bootstrap(simulated_genotype_matrix, simulated_variable_positions, 1000)

# Save each bootstrapped distribution as a gzipped csv.
results_prefix = './non_iua_models/ragsdale_2019/n_100/{0}/bootstraps/rep_id_{1}_'.format(f, rep_id)

np.savetxt(
    results_prefix+'ceu_abba.csv.gz',
    [bs_ceu[0]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_baba.csv.gz',
    [bs_ceu[1]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_bbaa.csv.gz',
    [bs_ceu[2]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_baaa.csv.gz',
    [bs_ceu[3]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_abaa.csv.gz',
    [bs_ceu[4]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_aaba.csv.gz',
    [bs_ceu[5]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_abba.csv.gz',
    [bs_chb[0]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_baba.csv.gz',
    [bs_chb[1]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_bbaa.csv.gz',
    [bs_chb[2]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_baaa.csv.gz',
    [bs_chb[3]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_abaa.csv.gz',
    [bs_chb[4]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_aaba.csv.gz',
    [bs_chb[5]],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)