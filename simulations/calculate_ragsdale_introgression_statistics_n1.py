### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = the admixture proportion ###


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

# Define Patterson's D function.
def pattersons_d(
    abba,
    baba,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA and BABA counts.
    ---------------------------------------------------------------------------
    OUTPUT: Patterson's D value.
    ###########################################################################
    """
    # Calculate Patterson's D.
    d_num = (abba - baba)
    d_den = (abba + baba)
    if (d_den != 0):
        d = (d_num / d_den)
    else:
        d = np.nan
    return d

# Define Danc function.
def danc(
    baaa,
    abaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide BAAA and ABAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: Danc value.
    ###########################################################################
    """
    # Calculate Danc.
    danc_num = (baaa - abaa)
    danc_den = (baaa + abaa)
    if (danc_den != 0):
        danc = (danc_num / danc_den)
    else:
        danc = np.nan
    return danc

# Define D+ function.
def dplus(
    abba,
    baba,
    baaa,
    abaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA, BABA, BAAA, and ABAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: D+ value.
    ###########################################################################
    """
    # Calculate D+.
    dplus_num = ((abba - baba) + (baaa - abaa))
    dplus_den = ((abba + baba) + (baaa + abaa))
    if (dplus_den != 0):
        dplus = (dplus_num / dplus_den)
    else:
        dplus = np.nan
    return dplus

# Define fhom function.
def fhom(
    abba,
    baba,
    aaba,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA, BABA, and AABA counts.
    ---------------------------------------------------------------------------
    OUTPUT: fhom value.
    ---------------------------------------------------------------------------
    NOTE: This equation only works for the case of one P3 sample.
    ###########################################################################
    """
    # Calculate fhom.
    fhom_num = (abba - baba)
    fhom_den = (abba + aaba)
    if (fhom_den != 0):
        fhom = (fhom_num / fhom_den)
    else:
        fhom = np.nan
    return fhom

# Define fanc function.
def fanc(
    baaa,
    abaa,
    bbaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide BAAA, ABAA, and BBAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: fanc value.
    ---------------------------------------------------------------------------
    NOTE: This equation only works for the case of one P3 sample.
    ###########################################################################
    """
    # Calculate fanc.
    fanc_num = (baaa - abaa)
    fanc_den = (baaa + bbaa)
    if (fanc_den != 0):
        fanc = (fanc_num / fanc_den)
    else:
        fanc = np.nan
    return fanc

# Define f+ function.
def fplus(
    abba,
    baba,
    aaba,
    baaa,
    abaa,
    bbaa,
):
    """
    ###########################################################################
    INPUT: Genome-wide ABBA, BABA, AABA, BAAA, ABAA, and BBAA counts.
    ---------------------------------------------------------------------------
    OUTPUT: f+ value.
    ---------------------------------------------------------------------------
    NOTE: This equation only works for the case of one P3 sample.
    ###########################################################################
    """
    # Calculate f+.
    fplus_num = ((abba - baba) + (baaa - abaa))
    fplus_den = ((abba + aaba) + (baaa + bbaa))
    if (fplus_den != 0):
        fplus = (fplus_num / fplus_den)
    else:
        fplus = np.nan
    return fplus


# Intialize arrays to store observed values.
# CEU.
ceu_obs_abba = np.array([])
ceu_obs_baba = np.array([])
ceu_obs_bbaa = np.array([])
ceu_obs_baaa = np.array([])
ceu_obs_abaa = np.array([])
ceu_obs_aaba = np.array([])
ceu_obs_d = np.array([])
ceu_obs_danc = np.array([])
ceu_obs_dplus = np.array([])
ceu_obs_fhom = np.array([])
ceu_obs_fanc = np.array([])
ceu_obs_fplus = np.array([])
# CHB.
chb_obs_abba = np.array([])
chb_obs_baba = np.array([])
chb_obs_bbaa = np.array([])
chb_obs_baaa = np.array([])
chb_obs_abaa = np.array([])
chb_obs_aaba = np.array([])
chb_obs_d = np.array([])
chb_obs_danc = np.array([])
chb_obs_dplus = np.array([])
chb_obs_fhom = np.array([])
chb_obs_fanc = np.array([])
chb_obs_fplus = np.array([])

# Parse comand-line arguments.
f = float(sys.argv[1])

# Calculate introgression statistics for all replicate simulations.
for rep_id in range(1, 101):
    # Load the genotype matrix.
    rep_genotype_matrix = np.loadtxt(
        './non_iua_models/ragsdale_2019/n_1/{0}/geno_mats/rep_id_{1}_geno_mat.csv.gz'.format(f, rep_id),
        dtype=int, delimiter=',',
    )
    # Calculate site patterns.
    ceu_rep_abba, ceu_rep_baba, ceu_rep_bbaa, ceu_rep_baaa, ceu_rep_abaa, ceu_rep_aaba = site_patterns(
        genotype_matrix=rep_genotype_matrix,
        p1_idx=[0],
        p2_idx=[1],
        p3_idx=[3],
    )
    # Calculate Patterson's D.
    ceu_rep_d = pattersons_d(
        abba=ceu_rep_abba,
        baba=ceu_rep_baba,
    )
    # Calculate Danc.
    ceu_rep_danc = danc(
        baaa=ceu_rep_baaa,
        abaa=ceu_rep_abaa,
    )
    # Calculate D+.
    ceu_rep_dplus = dplus(
        abba=ceu_rep_abba,
        baba=ceu_rep_baba,
        baaa=ceu_rep_baaa,
        abaa=ceu_rep_abaa,
    )
    # Calculate fhom.
    ceu_rep_fhom = fhom(
        abba=ceu_rep_abba,
        baba=ceu_rep_baba,
        aaba=ceu_rep_aaba,
    )
    # Calculate fanc.
    ceu_rep_fanc = fanc(
        baaa=ceu_rep_baaa,
        abaa=ceu_rep_abaa,
        bbaa=ceu_rep_bbaa,
    )
    # Calculate fplus.
    ceu_rep_fplus = fplus(
        abba=ceu_rep_abba,
        baba=ceu_rep_baba,
        aaba=ceu_rep_aaba,
        baaa=ceu_rep_baaa,
        abaa=ceu_rep_abaa,
        bbaa=ceu_rep_bbaa,
    )
    # Record the observed values for this ceu_replicate.
    ceu_obs_abba = np.append(ceu_obs_abba, ceu_rep_abba)
    ceu_obs_baba = np.append(ceu_obs_baba, ceu_rep_baba)
    ceu_obs_bbaa = np.append(ceu_obs_bbaa, ceu_rep_bbaa)
    ceu_obs_baaa = np.append(ceu_obs_baaa, ceu_rep_baaa)
    ceu_obs_abaa = np.append(ceu_obs_abaa, ceu_rep_abaa)
    ceu_obs_aaba = np.append(ceu_obs_aaba, ceu_rep_aaba)
    ceu_obs_d = np.append(ceu_obs_d, ceu_rep_d)
    ceu_obs_danc = np.append(ceu_obs_danc, ceu_rep_danc)
    ceu_obs_dplus = np.append(ceu_obs_dplus, ceu_rep_dplus)
    ceu_obs_fhom = np.append(ceu_obs_fhom, ceu_rep_fhom)
    ceu_obs_fanc = np.append(ceu_obs_fanc, ceu_rep_fanc)
    ceu_obs_fplus = np.append(ceu_obs_fplus, ceu_rep_fplus)
    # Calculate site patterns.
    chb_rep_abba, chb_rep_baba, chb_rep_bbaa, chb_rep_baaa, chb_rep_abaa, chb_rep_aaba = site_patterns(
        genotype_matrix=rep_genotype_matrix,
        p1_idx=[0],
        p2_idx=[2],
        p3_idx=[3],
    )
    # Calculate Patterson's D.
    chb_rep_d = pattersons_d(
        abba=chb_rep_abba,
        baba=chb_rep_baba,
    )
    # Calculate Danc.
    chb_rep_danc = danc(
        baaa=chb_rep_baaa,
        abaa=chb_rep_abaa,
    )
    # Calculate D+.
    chb_rep_dplus = dplus(
        abba=chb_rep_abba,
        baba=chb_rep_baba,
        baaa=chb_rep_baaa,
        abaa=chb_rep_abaa,
    )
    # Calculate fhom.
    chb_rep_fhom = fhom(
        abba=chb_rep_abba,
        baba=chb_rep_baba,
        aaba=chb_rep_aaba,
    )
    # Calculate fanc.
    chb_rep_fanc = fanc(
        baaa=chb_rep_baaa,
        abaa=chb_rep_abaa,
        bbaa=chb_rep_bbaa,
    )
    # Calculate fplus.
    chb_rep_fplus = fplus(
        abba=chb_rep_abba,
        baba=chb_rep_baba,
        aaba=chb_rep_aaba,
        baaa=chb_rep_baaa,
        abaa=chb_rep_abaa,
        bbaa=chb_rep_bbaa,
    )
    # Record the observed values for this chb_replicate.
    chb_obs_abba = np.append(chb_obs_abba, chb_rep_abba)
    chb_obs_baba = np.append(chb_obs_baba, chb_rep_baba)
    chb_obs_bbaa = np.append(chb_obs_bbaa, chb_rep_bbaa)
    chb_obs_baaa = np.append(chb_obs_baaa, chb_rep_baaa)
    chb_obs_abaa = np.append(chb_obs_abaa, chb_rep_abaa)
    chb_obs_aaba = np.append(chb_obs_aaba, chb_rep_aaba)
    chb_obs_d = np.append(chb_obs_d, chb_rep_d)
    chb_obs_danc = np.append(chb_obs_danc, chb_rep_danc)
    chb_obs_dplus = np.append(chb_obs_dplus, chb_rep_dplus)
    chb_obs_fhom = np.append(chb_obs_fhom, chb_rep_fhom)
    chb_obs_fanc = np.append(chb_obs_fanc, chb_rep_fanc)
    chb_obs_fplus = np.append(chb_obs_fplus, chb_rep_fplus)

# Save each observed result as a gzipped csv.
results_prefix = './non_iua_models/ragsdale_2019/n_1/{0}/obs_vals/'.format(f)

np.savetxt(
    results_prefix+'ceu_abba.csv.gz',
    [ceu_obs_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_baba.csv.gz',
    [ceu_obs_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_bbaa.csv.gz',
    [ceu_obs_bbaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_baaa.csv.gz',
    [ceu_obs_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_abaa.csv.gz',
    [ceu_obs_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_aaba.csv.gz',
    [ceu_obs_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_d.csv.gz',
    [ceu_obs_d],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_danc.csv.gz',
    [ceu_obs_danc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_dplus.csv.gz',
    [ceu_obs_dplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_fhom.csv.gz',
    [ceu_obs_fhom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_fanc.csv.gz',
    [ceu_obs_fanc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_fplus.csv.gz',
    [ceu_obs_fplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_abba.csv.gz',
    [chb_obs_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_baba.csv.gz',
    [chb_obs_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_bbaa.csv.gz',
    [chb_obs_bbaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_baaa.csv.gz',
    [chb_obs_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_abaa.csv.gz',
    [chb_obs_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_aaba.csv.gz',
    [chb_obs_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_d.csv.gz',
    [chb_obs_d],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_danc.csv.gz',
    [chb_obs_danc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_dplus.csv.gz',
    [chb_obs_dplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_fhom.csv.gz',
    [chb_obs_fhom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_fanc.csv.gz',
    [chb_obs_fanc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_fplus.csv.gz',
    [chb_obs_fplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)