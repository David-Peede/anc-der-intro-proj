### Dependencies ###
import multiprocessing
import numpy as np
import sys
### sys.argv[1] = the admixture proportion ###
### sys.argv[2] = replicate ID ###


# Define the bootstrapping class.
class Bootstrap:
    
    # Define static variables.
    simulated_genotype_matrix = np.loadtxt(
        './non_iua_models/ragsdale_2019/n_100/{0}/geno_mats/rep_id_{1}_geno_mat.csv.gz'.format(float(sys.argv[1]), int(sys.argv[2])),
        dtype=int, delimiter=',',
    )
    simulated_variable_positions = np.loadtxt(
        './non_iua_models/ragsdale_2019/n_100/{0}/var_pos/rep_id_{1}_var_pos.csv.gz'.format(float(sys.argv[1]), int(sys.argv[2])),
        delimiter=',',
    )
    
    # Define a derived allele frequency function.
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
        p1 = Bootstrap.derived_allele_freq(genotype_matrix, p1_idx)
        p2 = Bootstrap.derived_allele_freq(genotype_matrix, p2_idx)
        p3 = Bootstrap.derived_allele_freq(genotype_matrix, p3_idx)
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
    
    # Define a function to work on an instance of the window class.
    def do_bootstrapping(self):
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
        for i in range(100):
            # Intialize arrays to store bootstrapped values.
            ceu_results = np.zeros(6)
            chb_results = np.zeros(6)
            # For each bootstrap window (1000 x 100 kb).
            for j in range(1000):
                # Randomly generate a start and end position.
                start = np.random.randint(99_900_001)
                end = start + 100_000
                # Identify the variants that fall within the 100 kb window.
                variants = np.where(((start <= Bootstrap.simulated_variable_positions) & (Bootstrap.simulated_variable_positions <= end)))[0]
                # If there are variants to perform calculations on.
                if (variants.size > 0):
                    # Subset the 100 kb window from the genome-wide genotype matrix.
                    windowed_matrix = Bootstrap.simulated_genotype_matrix[variants, :]
                    # Calculate site pattern counts for the bootstrapped replicate and store the results.
                    # CEU.
                    ceu_abba, ceu_baba, ceu_baaa, ceu_abaa, _, _, _, _, ceu_bbaa, ceu_aaba = Bootstrap.site_patterns(
                        genotype_matrix=windowed_matrix,
                        p1_idx=np.arange(0, 100),
                        p2_idx=np.arange(100, 200),
                        p3_idx=np.arange(300, 301),
                    )
                    # CHB.
                    chb_abba, chb_baba, chb_baaa, chb_abaa, _, _, _, _, chb_bbaa, chb_aaba = Bootstrap.site_patterns(
                        genotype_matrix=windowed_matrix,
                        p1_idx=np.arange(0, 100),
                        p2_idx=np.arange(200, 300),
                        p3_idx=np.arange(300, 301),
                    )
                    # Append the results.
                    ceu_results += np.array([
                        ceu_abba, ceu_baba, ceu_bbaa,
                        ceu_baaa, ceu_abaa, ceu_aaba,
                    ])
                    chb_results += np.array([
                        chb_abba, chb_baba, chb_bbaa,
                        chb_baaa, chb_abaa, chb_aaba,
                    ])
                # Else...
                else:
                    # Continue to the next window.
                    continue
            # Append the results.
            ceu_abba_array = np.append(ceu_abba_array, ceu_results[0])
            ceu_baba_array = np.append(ceu_baba_array, ceu_results[1])
            ceu_bbaa_array = np.append(ceu_bbaa_array, ceu_results[2])
            ceu_baaa_array = np.append(ceu_baaa_array, ceu_results[3])
            ceu_abaa_array = np.append(ceu_abaa_array, ceu_results[4])
            ceu_aaba_array = np.append(ceu_aaba_array, ceu_results[5])
            chb_abba_array = np.append(chb_abba_array, chb_results[0])
            chb_baba_array = np.append(chb_baba_array, chb_results[1])
            chb_bbaa_array = np.append(chb_bbaa_array, chb_results[2])
            chb_baaa_array = np.append(chb_baaa_array, chb_results[3])
            chb_abaa_array = np.append(chb_abaa_array, chb_results[4])
            chb_aaba_array = np.append(chb_aaba_array, chb_results[5])
        # Aggregate results.
        ceu_bs_results = [ceu_abba_array, ceu_baba_array, ceu_bbaa_array, ceu_baaa_array, ceu_abaa_array, ceu_aaba_array]
        chb_bs_results = [chb_abba_array, chb_baba_array, chb_bbaa_array, chb_baaa_array, chb_abaa_array, chb_aaba_array]
        results_dicc = {'CEU': ceu_bs_results, 'CHB': chb_bs_results}
        return results_dicc
    
    # Define a function to do the work and return the bootstrapping results.
    def bootstrap_worker(bootstrap):
        return bootstrap.do_bootstrapping()

# Initialize an empty list to store instances of bootstrap.
bootstraps = []
# For 10 iterations...
for _ in range(10):
    # Create a bootstrap object and append the object list.
    bootstraps.append(Bootstrap())
# Pool as many process as cpus.
pool = multiprocessing.Pool(multiprocessing.cpu_count())
# Assign instances of bootstraps to pooled processes.
results = pool.map(Bootstrap.bootstrap_worker, bootstraps)
# Clean up from multiprocessing.
pool.close()
pool.join()
# Intialize a dictionaries to store results.
bs_ceu = {
    'abba': np.array([]), 'baba': np.array([]), 'bbaa': np.array([]),
    'baaa': np.array([]), 'abaa': np.array([]), 'aaba': np.array([]),
}
bs_chb = {
    'abba': np.array([]), 'baba': np.array([]), 'bbaa': np.array([]),
    'baaa': np.array([]), 'abaa': np.array([]), 'aaba': np.array([]),
}
# For all the result dictionaries...
for dicc in results:
    # Fill the dictionaries.
    bs_ceu['abba'] = np.append(bs_ceu['abba'], dicc['CEU'][0])
    bs_ceu['baba'] = np.append(bs_ceu['baba'], dicc['CEU'][1])
    bs_ceu['bbaa'] = np.append(bs_ceu['bbaa'], dicc['CEU'][2])
    bs_ceu['baaa'] = np.append(bs_ceu['baaa'], dicc['CEU'][3])
    bs_ceu['abaa'] = np.append(bs_ceu['abaa'], dicc['CEU'][4])
    bs_ceu['aaba'] = np.append(bs_ceu['aaba'], dicc['CEU'][5])
    bs_chb['abba'] = np.append(bs_chb['abba'], dicc['CHB'][0])
    bs_chb['baba'] = np.append(bs_chb['baba'], dicc['CHB'][1])
    bs_chb['bbaa'] = np.append(bs_chb['bbaa'], dicc['CHB'][2])
    bs_chb['baaa'] = np.append(bs_chb['baaa'], dicc['CHB'][3])
    bs_chb['abaa'] = np.append(bs_chb['abaa'], dicc['CHB'][4])
    bs_chb['aaba'] = np.append(bs_chb['aaba'], dicc['CHB'][5])

# Save each bootstrapped distribution as a gzipped csv.
results_prefix = './non_iua_models/ragsdale_2019/n_100/{0}/bootstraps/rep_id_{1}_'.format(float(sys.argv[1]), int(sys.argv[2]))

np.savetxt(
    results_prefix+'ceu_abba.csv.gz',
    [bs_ceu['abba']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_baba.csv.gz',
    [bs_ceu['baba']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_bbaa.csv.gz',
    [bs_ceu['bbaa']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_baaa.csv.gz',
    [bs_ceu['baaa']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_abaa.csv.gz',
    [bs_ceu['abaa']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'ceu_aaba.csv.gz',
    [bs_ceu['aaba']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_abba.csv.gz',
    [bs_chb['abba']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_baba.csv.gz',
    [bs_chb['baba']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_bbaa.csv.gz',
    [bs_chb['bbaa']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_baaa.csv.gz',
    [bs_chb['baaa']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_abaa.csv.gz',
    [bs_chb['abaa']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)

np.savetxt(
    results_prefix+'chb_aaba.csv.gz',
    [bs_chb['aaba']],
    fmt='%1.5f',
    delimiter=',',
    newline='\n',
)