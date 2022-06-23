import allel
import numcodecs
import numpy as np
import sys
import zarr


vcf_path = './zarr_arrays/tgp_aa_calls_alt_merged_filtered_chr{0}.vcf.gz'.format(str(sys.argv[1]))
zarr_path = './zarr_arrays/tgp_aa_calls_alt_merged_filtered_chr{0}.zarr'.format(str(sys.argv[1]))
allel.vcf_to_zarr(vcf_path, zarr_path, group=str(sys.argv[1]), fields='*', log=sys.stdout, overwrite=True)


callset = zarr.open_group(zarr_path, mode='r')
gt_zarr = callset[str(sys.argv[1])+'/calldata/GT']
print(gt_zarr.info)
gt = allel.GenotypeArray(gt_zarr)
print(gt)
print(gt.shape)
pos_sort = allel.SortedIndex(callset[str(sys.argv[1])+'/variants/POS'])
print(pos_sort)
print(pos_sort.shape)