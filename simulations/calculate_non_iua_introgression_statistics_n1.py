### Dependencies ###
import numpy as np
import sys
### sys.argv[1] = Neanderthal admixture proportion ###
### sys.argv[2] = Denisovan admixture proportion ###
### sys.argv[3] = non-iua model ###


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
# EUR.
eur_obs_abba = np.array([])
eur_obs_baba = np.array([])
eur_obs_bbaa = np.array([])
eur_obs_baaa = np.array([])
eur_obs_abaa = np.array([])
eur_obs_aaba = np.array([])
eur_obs_d = np.array([])
eur_obs_danc = np.array([])
eur_obs_dplus = np.array([])
eur_obs_fhom = np.array([])
eur_obs_fanc = np.array([])
eur_obs_fplus = np.array([])
# ASN.
asn_obs_abba = np.array([])
asn_obs_baba = np.array([])
asn_obs_bbaa = np.array([])
asn_obs_baaa = np.array([])
asn_obs_abaa = np.array([])
asn_obs_aaba = np.array([])
asn_obs_d = np.array([])
asn_obs_danc = np.array([])
asn_obs_dplus = np.array([])
asn_obs_fhom = np.array([])
asn_obs_fanc = np.array([])
asn_obs_fplus = np.array([])

# Parse comand-line arguments.
f_nea     = float(sys.argv[1])
f_den     = float(sys.argv[2])
non_iua   = str(sys.argv[3])

# Calculate introgression statistics for all replicate simulations.
for rep_id in range(1, 101):
    # Load the genotype matrix.
    rep_genotype_matrix = np.loadtxt(
        './non_iua_models/{0}/n_1/{1}_{2}/geno_mats/rep_id_{3}_geno_mat.csv.gz'.format(non_iua, f_nea, f_den, rep_id),
        dtype=int, delimiter=',',
    )
    # Calculate site patterns.
    eur_rep_abba, eur_rep_baba, eur_rep_bbaa, eur_rep_baaa, eur_rep_abaa, eur_rep_aaba = site_patterns(
        genotype_matrix=rep_genotype_matrix,
        p1_idx=[0],
        p2_idx=[1],
        p3_idx=[3],
    )
    # Calculate Patterson's D.
    eur_rep_d = pattersons_d(
        abba=eur_rep_abba,
        baba=eur_rep_baba,
    )
    # Calculate Danc.
    eur_rep_danc = danc(
        baaa=eur_rep_baaa,
        abaa=eur_rep_abaa,
    )
    # Calculate D+.
    eur_rep_dplus = dplus(
        abba=eur_rep_abba,
        baba=eur_rep_baba,
        baaa=eur_rep_baaa,
        abaa=eur_rep_abaa,
    )
    # Calculate fhom.
    eur_rep_fhom = fhom(
        abba=eur_rep_abba,
        baba=eur_rep_baba,
        aaba=eur_rep_aaba,
    )
    # Calculate fanc.
    eur_rep_fanc = fanc(
        baaa=eur_rep_baaa,
        abaa=eur_rep_abaa,
        bbaa=eur_rep_bbaa,
    )
    # Calculate fplus.
    eur_rep_fplus = fplus(
        abba=eur_rep_abba,
        baba=eur_rep_baba,
        aaba=eur_rep_aaba,
        baaa=eur_rep_baaa,
        abaa=eur_rep_abaa,
        bbaa=eur_rep_bbaa,
    )
    # Record the observed values for this eur_replicate.
    eur_obs_abba = np.append(eur_obs_abba, eur_rep_abba)
    eur_obs_baba = np.append(eur_obs_baba, eur_rep_baba)
    eur_obs_bbaa = np.append(eur_obs_bbaa, eur_rep_bbaa)
    eur_obs_baaa = np.append(eur_obs_baaa, eur_rep_baaa)
    eur_obs_abaa = np.append(eur_obs_abaa, eur_rep_abaa)
    eur_obs_aaba = np.append(eur_obs_aaba, eur_rep_aaba)
    eur_obs_d = np.append(eur_obs_d, eur_rep_d)
    eur_obs_danc = np.append(eur_obs_danc, eur_rep_danc)
    eur_obs_dplus = np.append(eur_obs_dplus, eur_rep_dplus)
    eur_obs_fhom = np.append(eur_obs_fhom, eur_rep_fhom)
    eur_obs_fanc = np.append(eur_obs_fanc, eur_rep_fanc)
    eur_obs_fplus = np.append(eur_obs_fplus, eur_rep_fplus)
    # Calculate site patterns.
    asn_rep_abba, asn_rep_baba, asn_rep_bbaa, asn_rep_baaa, asn_rep_abaa, asn_rep_aaba = site_patterns(
        genotype_matrix=rep_genotype_matrix,
        p1_idx=[0],
        p2_idx=[2],
        p3_idx=[3],
    )
    # Calculate Patterson's D.
    asn_rep_d = pattersons_d(
        abba=asn_rep_abba,
        baba=asn_rep_baba,
    )
    # Calculate Danc.
    asn_rep_danc = danc(
        baaa=asn_rep_baaa,
        abaa=asn_rep_abaa,
    )
    # Calculate D+.
    asn_rep_dplus = dplus(
        abba=asn_rep_abba,
        baba=asn_rep_baba,
        baaa=asn_rep_baaa,
        abaa=asn_rep_abaa,
    )
    # Calculate fhom.
    asn_rep_fhom = fhom(
        abba=asn_rep_abba,
        baba=asn_rep_baba,
        aaba=asn_rep_aaba,
    )
    # Calculate fanc.
    asn_rep_fanc = fanc(
        baaa=asn_rep_baaa,
        abaa=asn_rep_abaa,
        bbaa=asn_rep_bbaa,
    )
    # Calculate fplus.
    asn_rep_fplus = fplus(
        abba=asn_rep_abba,
        baba=asn_rep_baba,
        aaba=asn_rep_aaba,
        baaa=asn_rep_baaa,
        abaa=asn_rep_abaa,
        bbaa=asn_rep_bbaa,
    )
    # Record the observed values for this asn_replicate.
    asn_obs_abba = np.append(asn_obs_abba, asn_rep_abba)
    asn_obs_baba = np.append(asn_obs_baba, asn_rep_baba)
    asn_obs_bbaa = np.append(asn_obs_bbaa, asn_rep_bbaa)
    asn_obs_baaa = np.append(asn_obs_baaa, asn_rep_baaa)
    asn_obs_abaa = np.append(asn_obs_abaa, asn_rep_abaa)
    asn_obs_aaba = np.append(asn_obs_aaba, asn_rep_aaba)
    asn_obs_d = np.append(asn_obs_d, asn_rep_d)
    asn_obs_danc = np.append(asn_obs_danc, asn_rep_danc)
    asn_obs_dplus = np.append(asn_obs_dplus, asn_rep_dplus)
    asn_obs_fhom = np.append(asn_obs_fhom, asn_rep_fhom)
    asn_obs_fanc = np.append(asn_obs_fanc, asn_rep_fanc)
    asn_obs_fplus = np.append(asn_obs_fplus, asn_rep_fplus)

# Save each observed result as a gzipped csv.
results_prefix = './non_iua_models/{0}/n_1/{1}_{2}/obs_vals/'.format(non_iua, f_nea, f_den)

np.savetxt(
    results_prefix+'eur_abba.csv.gz',
    [eur_obs_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_baba.csv.gz',
    [eur_obs_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_bbaa.csv.gz',
    [eur_obs_bbaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_baaa.csv.gz',
    [eur_obs_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_abaa.csv.gz',
    [eur_obs_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_aaba.csv.gz',
    [eur_obs_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_d.csv.gz',
    [eur_obs_d],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_danc.csv.gz',
    [eur_obs_danc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_dplus.csv.gz',
    [eur_obs_dplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_fhom.csv.gz',
    [eur_obs_fhom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_fanc.csv.gz',
    [eur_obs_fanc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'eur_fplus.csv.gz',
    [eur_obs_fplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_abba.csv.gz',
    [asn_obs_abba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_baba.csv.gz',
    [asn_obs_baba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_bbaa.csv.gz',
    [asn_obs_bbaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_baaa.csv.gz',
    [asn_obs_baaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_abaa.csv.gz',
    [asn_obs_abaa],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_aaba.csv.gz',
    [asn_obs_aaba],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_d.csv.gz',
    [asn_obs_d],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_danc.csv.gz',
    [asn_obs_danc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_dplus.csv.gz',
    [asn_obs_dplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_fhom.csv.gz',
    [asn_obs_fhom],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_fanc.csv.gz',
    [asn_obs_fanc],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'asn_fplus.csv.gz',
    [asn_obs_fplus],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)