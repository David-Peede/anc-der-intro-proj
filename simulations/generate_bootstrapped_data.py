### Dependencies ###
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

# Define dXY function.
def dxy_freq(
    genotype_matrix,
    px_idx,
    py_idx,
    sequence_length=100_000_000,
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genotypes from msprime.
        px_idx: PX lineage.
        py_idx: PY lineage.
        sequence_length: Length of the genome simulated.
    ---------------------------------------------------------------------------
    OUTPUT: Mean pairwise sequence divergence between PX and PY across all
            sites.
    ###########################################################################
    """
    # Calculate derived allele frequencies per taxa.
    px = derived_allele_freq(genotype_matrix, px_idx)
    py = derived_allele_freq(genotype_matrix, py_idx)
    # Calculate sequence divergence.
    pxy = ((px) * (1 - py))
    pyx = ((py) * (1 - px))
    psum = (pxy + pyx)
    ptot = psum.sum()
    # Divide by the sequence length for average pairwise sequence divergence.
    dxy = (ptot / sequence_length)
    return dxy

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

# Define genetic distance function.
def genetic_distances(
    genotype_matrix,
    p1_idx=[0],
    p2_idx=[1],
    p3_idx=[2],
    sequence_length=100_000_000,
):
    """
    ###########################################################################
    INPUT
        genotype_matrix: Simulated genotypes from msprime.
        p1_idx: P1 lineage.
        p2_idx: P2 lineage.
        p3_idx: P3 lineage.
        sequence_length: Length of the genome simulated.
    ---------------------------------------------------------------------------
    OUTPUT: Mean pairwise sequence divergence between all combinations of taxa.
    ###########################################################################
    """
    # Calculate the genetic distance for all combinations of taxa.
    d12 = dxy_freq(
        genotype_matrix,
        p1_idx,
        p2_idx,
        sequence_length,
    )
    d13 = dxy_freq(
        genotype_matrix,
        p1_idx,
        p3_idx,
        sequence_length,
    )
    d23 = dxy_freq(
        genotype_matrix,
        p2_idx,
        p3_idx,
        sequence_length,
    )
    return d12, d13, d23

# Define bootstrapping function.
def genotype_matrix_bootstrap(
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
    abba_array = np.array([])
    baba_array = np.array([])
    bbaa_array = np.array([])
    baaa_array = np.array([])
    abaa_array = np.array([])
    aaba_array = np.array([])
    d12_array = np.array([])
    d13_array = np.array([])
    d23_array = np.array([])
    # For each bootstrap replicate.
    for i in range(n_replicates):
        # Intialize an empty array to store the bootstrapped genome.
        bootstrapped_genome = np.empty((0, 3))
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
        abba, baba, bbaa, baaa, abaa, aaba = site_patterns(
            genotype_matrix=bootstrapped_genome,
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
        # Calculate genetic distances for the bootstrapped replicate and store the results.
        d12, d13, d23 = genetic_distances(
            genotype_matrix=bootstrapped_genome,
            p1_idx=[0],
            p2_idx=[1],
            p3_idx=[2],
            sequence_length=100_000_000,
        )
        d12_array = np.append(d12_array, d12)
        d13_array = np.append(d13_array, d13)
        d23_array = np.append(d23_array, d23)
    return abba_array, baba_array, bbaa_array, baaa_array, abaa_array, aaba_array, d12_array, d13_array, d23_array


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

# Perform bootstrapping.
bs_abba, bs_baba, bs_bbaa, bs_baaa, bs_abaa, bs_aaba, bs_d12, bs_d13, bs_d23 = genotype_matrix_bootstrap(simulated_genotype_matrix, simulated_variable_positions, 1000)

# Save each bootstrapped distribution as a gzipped csv.
results_prefix = './sim_outputs/{0}/bootstraps/rep_id_{1}_'.format(f, rep_id)

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
    results_prefix+'bbaa.csv.gz',
    [bs_bbaa],
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
    results_prefix+'aaba.csv.gz',
    [bs_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'d12.csv.gz',
    [bs_d12],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'d13.csv.gz',
    [bs_d13],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'d23.csv.gz',
    [bs_d23],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)