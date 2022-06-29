### Dependencies ###
import gzip
import sys
### sys.argv[1] = unfiltered vcf ###


# Define function to filter the tgp joint vcf file.
def filter_tgp_chimp_nean_vcf(vcf):
    """
    ###########################################################################
    INPUT: A gzipped unfiltered VCF file where the last index is Altai
           Neanderthal and the second to last index is the chimp.
    ---------------------------------------------------------------------------
    OUTPUT: Filtered VCF file to standard out.
    ###########################################################################
    """
    # Intialize the filtered vcf file to be outputted to stdout.
    new_vcf = sys.stdout
    # Open the original vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line is a header line write it to the new vcf.
            if line.startswith('##') or line.startswith('#'):
                new_vcf.write(line)
            # Else split the line by tabs.
            else:
                spline = line.split()
                # Determine the length of the refernce and alternative fields.
                ref_len = len(spline[3])
                alt_len = len(spline[4])
                # If the length of refernce and alternative fields is not
                # exactly 2, ie if the site is not mono- or biallelic,
                # continue to the next line.
                if ((ref_len + alt_len) != 2):
                    continue
                # Else grab the chimp genotype info for this site.
                else:
                    chimp_gt = spline[-2]
                    # If there isn't genotype information for the chimp
                    # continue to the next line.
                    if (chimp_gt == './.'):
                        continue
                    # Else grab the neanderthal genotype info for this site.
                    else:
                        nean_gt = spline[-1]
                        # If the genotype information for the neanderthal is
                        # missing continue to the next line.
                        if (nean_gt[0] == '.'):
                            continue
                        # Else grab the genotype quality score for the
                        # neanderthal.
                        else:
                            nean_gq = int(nean_gt.split(':')[-1])
                            # If the genotype quality score is less than 40
                            # continue to the next line.
                            if (nean_gq < 40):
                                continue
                            # Else grab the genotype info for the first
                            # individual in the tgp.
                            else:
                                tgp_gt = spline[9]
                                # If the tgp has an allele call then write
                                # the line to the new vcf file.
                                if (tgp_gt[0] != '.'):
                                    new_vcf.write(line)
                                # Else make all of the tgp individuals
                                # homozygous for the reference allele since
                                # the tgp data is imputed.
                                else:
                                    new_spline = [geno.replace('./.:.:.:.:.:.:.:.', '0|0:.:.:.:.:.:.:.') for geno in spline]
                                    new_vcf.write('\t'.join(new_spline)+'\n')
    return


filter_tgp_chimp_nean_vcf(vcf=str(sys.argv[1]))