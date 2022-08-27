# Import packages.
import gzip
import numpy as np


# Define a function to record the positions in a vcf file.
def canid_vcf_positions():
    """
    ###########################################################################
    INPUT: N/A
    ---------------------------------------------------------------------------
    OUTPUT: Gzipped csv containg all the positions present in the vcf file.
    ###########################################################################
    """
    # For every chromosome...
    for chrom in range(1, 39):
        # Intialize an empty array to store positions.
        pos_array = np.array([])
        # Grab the vcf file.
        vcf = './canid_merged_filtered_chr{0}.vcf.gz'.format(chrom)
        # Using the gzip package open the vcf file.
        with gzip.open(vcf, 'rt') as data:
            # Loop through every line in the vcf file.
            for line in data:
                # If the current line is a part of the meta info...
                if line.startswith('#'):
                    # Continue to the next line in the vcf file.
                    continue
                # Else...
                else:
                    # Split the line.
                    spline = line.split()
                    # Grab the current position.
                    pos = int(spline[1])
                    # Append the positions array.
                    pos_array = np.append(pos_array, pos)
        # Save the positions array.
        np.savetxt(
            './canid_positions_chr{0}.csv.gz'.format(chrom),
            [pos_array],
            fmt='%d',
            delimiter=',',
            )
    return


# Determine the positions that passed QC for each chromosome.
canid_vcf_positions()