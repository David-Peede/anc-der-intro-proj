### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = Neanderthal admixture proportion ###
### sys.argv[2] = Denisovan admixture proportion ###
### sys.argv[3] = non-iua model ###
### sys.argv[4] = replicate ID ###
### sys.argv[5] = sample size ###


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
f_nea     = float(sys.argv[1])
f_den     = float(sys.argv[2])
non_iua   = str(sys.argv[3])
rep_id    = int(sys.argv[4])

# Load in the replicate data.
simulated_genotype_matrix = np.loadtxt(
    './non_iua_models/{0}/n_100/{1}_{2}/geno_mats/rep_id_{3}_geno_mat.csv.gz'.format(non_iua, f_nea, f_den, rep_id),
    dtype=int, delimiter=',',
)

# Calculate site patterns.
eur_abba, eur_baba, eur_baaa, eur_abaa, eur_abba_hom, eur_baba_hom, eur_baaa_hom, eur_abaa_hom, eur_bbaa, eur_aaba = site_patterns(
    genotype_matrix=simulated_genotype_matrix,
    p1_idx=np.arange(0, 100),
    p2_idx=np.arange(100, 200),
    p3_idx=np.arange(300, 301),
)
# Calculate detection metrics.
eur_d, eur_danc, eur_dplus = detection_metrics(
    abba=eur_abba,
    baba=eur_baba,
    baaa=eur_baaa,
    abaa=eur_abaa,
)
# Calculate quantification metrics.
eur_fhom, eur_fanc, eur_fplus = quantification_metrics(
    abba_num=eur_abba,
    baba_num=eur_baba,
    baaa_num=eur_baaa,
    abaa_num=eur_abaa,
    abba_den=eur_abba_hom,
    baba_den=eur_baba_hom,
    baaa_den=eur_baaa_hom,
    abaa_den=eur_abaa_hom,
)
# Consolidate the results.
eur_site_patterns = np.array([eur_abba, eur_baba, eur_bbaa, eur_baaa, eur_abaa, eur_aaba])
eur_detection = np.array([eur_d, eur_danc, eur_dplus])
eur_quantification = np.array([eur_fhom, eur_fanc, eur_fplus])
# Calculate site patterns.
asn_abba, asn_baba, asn_baaa, asn_abaa, asn_abba_hom, asn_baba_hom, asn_baaa_hom, asn_abaa_hom, asn_bbaa, asn_aaba = site_patterns(
    genotype_matrix=simulated_genotype_matrix,
    p1_idx=np.arange(0, 100),
    p2_idx=np.arange(200, 300),
    p3_idx=np.arange(300, 301),
)
# Calculate detection metrics.
asn_d, asn_danc, asn_dplus = detection_metrics(
    abba=asn_abba,
    baba=asn_baba,
    baaa=asn_baaa,
    abaa=asn_abaa,
)
# Calculate quantification metrics.
asn_fhom, asn_fanc, asn_fplus = quantification_metrics(
    abba_num=asn_abba,
    baba_num=asn_baba,
    baaa_num=asn_baaa,
    abaa_num=asn_abaa,
    abba_den=asn_abba_hom,
    baba_den=asn_baba_hom,
    baaa_den=asn_baaa_hom,
    abaa_den=asn_abaa_hom,
)
# Consolidate the results.
asn_site_patterns = np.array([asn_abba, asn_baba, asn_bbaa, asn_baaa, asn_abaa, asn_aaba])
asn_detection = np.array([asn_d, asn_danc, asn_dplus])
asn_quantification = np.array([asn_fhom, asn_fanc, asn_fplus])

# Save each observed result as a gzipped csv.
results_prefix = './non_iua_models/{0}/n_100/{1}_{2}/obs_vals/rep_id_{3}_'.format(non_iua, f_nea, f_den, rep_id)

np.savetxt(
    results_prefix+'eur_site_patterns.csv.gz',
    [eur_site_patterns],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_detection_metrics.csv.gz',
    [eur_detection],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_quantification_metrics.csv.gz',
    [eur_quantification],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_site_patterns.csv.gz',
    [asn_site_patterns],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_detection_metrics.csv.gz',
    [asn_detection],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_quantification_metrics.csv.gz',
    [asn_quantification],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
