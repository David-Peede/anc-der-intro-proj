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
    './non_iua_models/ragsdale_2019/n_100/{0}/geno_mats/rep_id_{1}_geno_mat.csv.gz'.format(f, rep_id),
    dtype=int, delimiter=',',
)

# Calculate site patterns.
ceu_abba, ceu_baba, ceu_baaa, ceu_abaa, ceu_abba_hom, ceu_baba_hom, ceu_baaa_hom, ceu_abaa_hom, ceu_bbaa, ceu_aaba = site_patterns(
    genotype_matrix=simulated_genotype_matrix,
    p1_idx=np.arange(0, 100),
    p2_idx=np.arange(100, 200),
    p3_idx=np.arange(300, 301),
)
# Calculate detection metrics.
ceu_d, ceu_danc, ceu_dplus = detection_metrics(
    abba=ceu_abba,
    baba=ceu_baba,
    baaa=ceu_baaa,
    abaa=ceu_abaa,
)
# Calculate quantification metrics.
ceu_fhom, ceu_fanc, ceu_fplus = quantification_metrics(
    abba_num=ceu_abba,
    baba_num=ceu_baba,
    baaa_num=ceu_baaa,
    abaa_num=ceu_abaa,
    abba_den=ceu_abba_hom,
    baba_den=ceu_baba_hom,
    baaa_den=ceu_baaa_hom,
    abaa_den=ceu_abaa_hom,
)
# Consolidate the results.
ceu_site_patterns = np.array([ceu_abba, ceu_baba, ceu_bbaa, ceu_baaa, ceu_abaa, ceu_aaba])
ceu_detection = np.array([ceu_d, ceu_danc, ceu_dplus])
ceu_quantification = np.array([ceu_fhom, ceu_fanc, ceu_fplus])
# Calculate site patterns.
chb_abba, chb_baba, chb_baaa, chb_abaa, chb_abba_hom, chb_baba_hom, chb_baaa_hom, chb_abaa_hom, chb_bbaa, chb_aaba = site_patterns(
    genotype_matrix=simulated_genotype_matrix,
    p1_idx=np.arange(0, 100),
    p2_idx=np.arange(200, 300),
    p3_idx=np.arange(300, 301),
)
# Calculate detection metrics.
chb_d, chb_danc, chb_dplus = detection_metrics(
    abba=chb_abba,
    baba=chb_baba,
    baaa=chb_baaa,
    abaa=chb_abaa,
)
# Calculate quantification metrics.
chb_fhom, chb_fanc, chb_fplus = quantification_metrics(
    abba_num=chb_abba,
    baba_num=chb_baba,
    baaa_num=chb_baaa,
    abaa_num=chb_abaa,
    abba_den=chb_abba_hom,
    baba_den=chb_baba_hom,
    baaa_den=chb_baaa_hom,
    abaa_den=chb_abaa_hom,
)
# Consolidate the results.
chb_site_patterns = np.array([chb_abba, chb_baba, chb_bbaa, chb_baaa, chb_abaa, chb_aaba])
chb_detection = np.array([chb_d, chb_danc, chb_dplus])
chb_quantification = np.array([chb_fhom, chb_fanc, chb_fplus])

# Save each observed result as a gzipped csv.
results_prefix = './non_iua_models/ragsdale_2019/n_100/{0}/obs_vals/rep_id_{1}_'.format(f, rep_id)

np.savetxt(
    results_prefix+'ceu_site_patterns.csv.gz',
    [ceu_site_patterns],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_detection_metrics.csv.gz',
    [ceu_detection],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_quantification_metrics.csv.gz',
    [ceu_quantification],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_site_patterns.csv.gz',
    [chb_site_patterns],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_detection_metrics.csv.gz',
    [chb_detection],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_quantification_metrics.csv.gz',
    [chb_quantification],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)
