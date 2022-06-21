import gzip
import re
import sys

#for i in range(1,23):
#    chr_list.append(str(i))
chr_list = [int(sys.argv[1])]

for chromosome in chr_list:
    infile='tgp_biallelic_snps_chr{0}.vcf.gz'.format(chromosome)
    outfile='tgp_biallelic_snps_chr{0}_aa_calls.vcf'.format(chromosome)
    g = open(outfile, 'w')
    with gzip.open(infile, 'rt') as f:
        for line in f:
            if '##' in line:
                g.write(line)
            elif '#' in line:
                spline = line.split()
                spline.append('Ancestor')
                outLine = '\t'.join(spline)
                g.write(outLine)
            else:
                spline = line.split()
                if 'AA' in spline[7] and 'VT=SNP' in spline[7]:
                    anc_allele_search = re.search('^.+;AA=(\S+)\|\|\|.+', spline[7])
                    anc_allele = anc_allele_search.group(1)
                    if spline[3] == anc_allele.upper():
                        anc_allele = '0|0'
                    elif spline[4] == anc_allele.upper():
                        anc_allele = '1|1'
                    else:
                        anc_allele = './.'
                else:
                    anc_allele = './.'
                spline.append(anc_allele)
                outLine = '\n'+'\t'.join(spline)
                g.write(outLine)
g.close()