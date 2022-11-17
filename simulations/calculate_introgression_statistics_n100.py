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

# Define a function to calculate site patterns.
def site_patterns(
    genotype_matrix,
    p1_idx,
    p2_idx,
    p3_idx,
):
    """
    ###########################################################################
    INPUT
        p1: Derived alle frequency array for the P1 lineage.
        p2: Derived alle frequency array for the P2 lineage.
        p3: Derived alle frequency array for the P3 lineage.
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
    baaa = ((p1) * (1 - p2) * (1 - p3)).sum()
    abaa = ((1 - p1) * (p2) * (1 - p3)).sum()
    bbaa = ((p1) * (p2) * (1 - p3)).sum()
    aaba = ((1 - p1) * (1 - p2) * (p3)).sum()
    abba_hom = ((1 - p1) * (p3) * (p3)).sum()
    baba_hom = ((p1) * (1 - p3) * (p3)).sum()
    baaa_hom = ((p1) * (1 - p3) * (1 - p3)).sum()
    abaa_hom = ((1 - p1) * (p3) * (1 - p3)).sum()
    return abba, baba, baaa, abaa, abba_hom, baba_hom, baaa_hom, abaa_hom, bbaa, aaba

# Define a function to calculate detection metrics.
def detection_metrics(
    abba,
    baba,
    baaa,
    abaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA, BABA, BAAA, and ABAA counts from frequencies.
    ---------------------------------------------------------------------------
    OUTPUT: Patterson's D, Danc, and D+ values.
    ###########################################################################
    """
    # Calculate Patterson's D.
    d = ((abba - baba) / (abba + baba))
    # Calculate Danc.
    danc = ((baaa - abaa) / (baaa + abaa))
    # Calculate D+.
    dplus = (((abba - baba) + (baaa - abaa)) / ((abba + baba) + (baaa + abaa)))
    return d, danc, dplus

# Define a function to calculate quantification metrics.
def quantification_metrics(
    abba_num,
    baba_num,
    baaa_num,
    abaa_num,
    abba_den,
    baba_den,
    baaa_den,
    abaa_den,
):
    """
    ###########################################################################
    INPUT
        {site pattern}_num: Site pattern counts from ((P1, P2), P3).
        {site pattern}_den: Site pattern counts from ((P1, P3), P3).
    ---------------------------------------------------------------------------
    OUTPUT: fhom, fanc, and f+ values.
    ###########################################################################
    """
    # Calculate fhom.
    fhom = ((abba_num - baba_num) / (abba_den - baba_den))
    # Calculate fanc.
    fanc = ((baaa_num - abaa_num) / (baaa_den - abaa_den))
    # Calculate f+.
    fplus = (((abba_num - baba_num) + (baaa_num - abaa_num)) / ((abba_den - baba_den) + (baaa_den - abaa_den)))
    return fhom, fanc, fplus

# Parse comand-line arguments.
f         = float(sys.argv[1])
rep_id    = int(sys.argv[2])

# Load in the replicate data.
simulated_genotype_matrix = np.loadtxt(
    './sim_outputs/n_100/{0}/geno_mats/rep_id_{1}_geno_mat.csv.gz'.format(f, rep_id),
    dtype=int, delimiter=',',
)

# Calculate site patterns.
obs_abba, obs_baba, obs_baaa, obs_abaa, obs_abba_hom, obs_baba_hom, obs_baaa_hom, obs_abaa_hom, obs_bbaa, obs_aaba = site_patterns(
    genotype_matrix=simulated_genotype_matrix,
    p1_idx=np.arange(0, 100),
    p2_idx=np.arange(100, 200),
    p3_idx=np.arange(200, 201),
)
# Calculate detection metrics.
obs_d, obs_danc, obs_dplus = detection_metrics(
    abba=obs_abba,
    baba=obs_baba,
    baaa=obs_baaa,
    abaa=obs_abaa,
)
# Calculate quantification metrics.
obs_fhom, obs_fanc, obs_fplus = quantification_metrics(
    abba_num=obs_abba,
    baba_num=obs_baba,
    baaa_num=obs_baaa,
    abaa_num=obs_abaa,
    abba_den=obs_abba_hom,
    baba_den=obs_baba_hom,
    baaa_den=obs_baaa_hom,
    abaa_den=obs_abaa_hom,
)
# Consolidate the results.
obs_site_patterns = np.array([obs_abba, obs_baba, obs_bbaa, obs_baaa, obs_abaa, obs_aaba])
obs_detection = np.array([obs_d, obs_danc, obs_dplus])
obs_quantification = np.array([obs_fhom, obs_fanc, obs_fplus])

# Save each observed result as a gzipped csv.
results_prefix = './sim_outputs/n_100/{0}/obs_vals/rep_id_{1}_'.format(f, rep_id)

np.savetxt(
    results_prefix+'obs_site_patterns.csv.gz',
    [obs_site_patterns],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'obs_detection_metrics.csv.gz',
    [obs_detection],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'obs_quantification_metrics.csv.gz',
    [obs_quantification],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
