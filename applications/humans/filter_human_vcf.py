### Dependencies ###
import gzip
import sys
### sys.argv[1] = unfiltered vcf ###
### sys.argv[2] = chromosome ###

# Define function to filter the tgp joint vcf file.
def filter_tgp_archaic_vcf(vcf, chrom):
    """
    ###########################################################################
    INPUT
        vcf: A gzipped unfiltered VCF file where the last index is the archaic.
        chrom: Chromosome for the VCF.
    ---------------------------------------------------------------------------
    OUTPUT: Filtered VCF file to standard out.
    ###########################################################################
    """
    # Load the fasta reference for the ancestral allele calls
    anc_fa = './epo_calls/homo_sapiens_ancestor_{0}.fa'.format(chrom)
    # Intialize an empty string for storing the data.
    anc_seq = ''
    # Open the fasta file.
    with open(anc_fa) as data:
        # For every line in data...
        for line in data:
            # If the line is a header...
            if line.startswith('>'):
                # Continue to the next line.
                continue
            # Else...
            else:
                # Append the line to the sequence string.
                anc_seq += line.strip().replace(' ','')
    # Intialize the filtered vcf file to be outputted to stdout.
    new_vcf = sys.stdout
    # Open the original vcf file.
    with gzip.open(vcf, 'rt') as data:
        # Iterate through every line in the original vcf file.
        for line in data:
            # If the line is a meta info line...
            if line.startswith('##'):
                # Write it to the new vcf.
                new_vcf.write(line)
            # Else-if the line is the header line...
            elif line.startswith('#'):
                # Split the line by tabs.
                spline = line.split()
                # Append the ancestor column.
                spline.append('Ancestor')
                # Write it to the new vcf.
                new_vcf.write('\t'.join(spline)+'\n')
            # Else...
            else:
                # Split the line by tabs.
                spline = line.split()
                # Grab the refernce and alternative alleles.
                ref = spline[3]
                alt = spline[4]
                # Grab the position.
                pos = int(spline[1])
                # Determine the ancestral allele.
                anc_allele = anc_seq[pos-1]
                # Grab the archaic site info.
                arc_site = spline[-1]
                # Grab the genotype info for the first individual in the tgp.
                tgp_gt = spline[9]
                # If the length of refernce and alternative fields is not
                # exactly 2, ie if the site is not mono- or biallelic...
                if ((len(ref) + len(alt)) != 2):
                    # Continue to the next line.
                    continue
                # Else-if the archaic site is missing...
                elif (arc_site[0] == '.'):
                    # Continue to the next line.
                    continue
                # Else-if the archaic has a GQ lower than 40..
                elif (int(arc_site.split(':')[-1]) < 40):
                    # Continue to the next line.
                    continue
                # Else-if the ancestral allele is not known.
                elif ((anc_allele == '.') | (anc_allele == '-') | (anc_allele == 'N')):
                    # Continue to the next line.
                    continue
                # Else-if the ancestral allele does not match the refernce or alternate
                # allele and there is an allele call for alternate allele...
                elif ((anc_allele.upper() != ref) & (anc_allele.upper() != alt) & (alt != '.')):
                    # Continue to the next line.
                    continue
                # Else...
                else:
                    # If the ancestral allele is the refernce allele...
                    if (anc_allele.upper() == ref):
                        # Set ancestral genotype to be homozygous refernce.
                        anc_gt = '0|0:.:.:.:.:.:.:.'
                        # If the tgp has an allele call...
                        if (tgp_gt[0] != '.'):
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(spline)+'\n')
                        # Else...
                        else:
                            # Append the info field.
                            spline[7] = spline[7]+';AA={0}|||'.format(anc_allele)
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Make all of the tgp individuals homozygous for the reference allele since
                            # the tgp data is imputed.
                            new_spline = [geno.replace('./.:.:.:.:.:.:.:.', '0|0:.:.:.:.:.:.:.') for geno in spline]
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(new_spline)+'\n')
                    # Else-if the ancestral allele is the alternative allele...
                    elif (anc_allele.upper() == alt):
                        # Set ancestral genotype to be homozygous alternative.
                        anc_gt = '1|1:.:.:.:.:.:.:.'
                        # If the tgp has an allele call...
                        if (tgp_gt[0] != '.'):
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(spline)+'\n')
                        # Else...
                        else:
                            # Append the info field.
                            spline[7] = spline[7]+';AA={0}|||'.format(anc_allele)
                            # Append the Ancestor column.
                            spline.append(anc_gt)
                            # Make all of the tgp individuals homozygous for the reference allele since
                            # the tgp data is imputed.
                            new_spline = [geno.replace('./.:.:.:.:.:.:.:.', '0|0:.:.:.:.:.:.:.') for geno in spline]
                            # Write it to the new vcf.
                            new_vcf.write('\t'.join(new_spline)+'\n')
                    # Else...
                    else:
                        # Set ancestral genotype to be homozygous alternative.
                        anc_gt = '1|1:.:.:.:.:.:.:.'
                        # Set the alternative allele as the ancestral allele.
                        spline[4] = anc_allele.upper()
                        # Append the info field.
                        spline[7] = spline[7]+';AA={0}|||'.format(anc_allele)
                        # Append the Ancestor column.
                        spline.append(anc_gt)
                        # Make all of the tgp individuals homozygous for the reference allele since
                        # the tgp data is imputed.
                        new_spline = [geno.replace('./.:.:.:.:.:.:.:.', '0|0:.:.:.:.:.:.:.') for geno in spline]
                        # Write it to the new vcf.
                        new_vcf.write('\t'.join(new_spline)+'\n') 
    return


filter_tgp_archaic_vcf(vcf=str(sys.argv[1]), chrom=str(sys.argv[2]))