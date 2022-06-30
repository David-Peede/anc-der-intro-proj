### Dependencies ###
import allel
import numcodecs
import numpy as np
import sys
import zarr
### sys.argv[1] = chromosome number ###


# Define file paths.
vcf_path = './zarr_arrays/tgp_chimp_altai_merged_filtered_chr{0}.vcf.gz'.format(str(sys.argv[1]))
zarr_path = './zarr_arrays/tgp_chimp_altai_merged_filtered_chr{0}.zarr'.format(str(sys.argv[1]))

# Convert the vcf file to a zarr array.
allel.vcf_to_zarr(vcf_path, zarr_path, group=str(sys.argv[1]), fields='*', log=sys.stdout, overwrite=True)

# Open the zarr array.
callset = zarr.open_group(zarr_path, mode='r')
# Extract the genotype information.
gt_zarr = callset[str(sys.argv[1])+'/calldata/GT']
# Print the genotype information.
print(gt_zarr.info)
# Convert the gentoype information to a genotype array.
gt = allel.GenotypeArray(gt_zarr)
# Print the genotype array.
print(gt)
print(gt.shape)
# Extract the positions for genotypes.
pos_sort = allel.SortedIndex(callset[str(sys.argv[1])+'/variants/POS'])
# Print the position information.
print(pos_sort)
print(pos_sort.shape)
# Check if there is any missing genotype information and if so how much.
gt_missing = gt.is_missing()
print(np.count_nonzero(gt_missing))