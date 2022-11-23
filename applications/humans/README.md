# Human Application

This directory contains all the code to replicate the human analyses from Peede et al. 202X, which I will now outline.

## Download

All the data used is publicly available and can be downloaded from the following locations:

- [Phase 3 Release of the 1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)
- [Altai Neanderthal](http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/)
- [EPO Ancestral Allele Calls](http://ftp.ensembl.org/pub/release-74/fasta/ancestral_alleles/)

## Merge

I used `bcftools v1.13` to merge the TGP and Altai Neanderthal data.

```bash
# Merge the TGP and Altai Neanderthal VCF files.
for CHR in {1..22}; do
bcftools merge -Oz -o tgp_altai_merged_raw_chr${CHR}.vcf.gz ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz chr${CHR}_mq25_mapab100.vcf.gz
done
```

## Filter

I used `filter_human_vcf.py` to create an all sites VCF file.

```bash
# Filter the merged VCF files.
for CHR in {1..22}; do
python filter_human_vcf.py tgp_altai_merged_raw_chr${CHR}.vcf.gz ${CHR} | bgzip > tgp_altai_merged_filtered_all_sites_chr${CHR}.vcf.gz
done
```

To save space you can optionally remove mono-allelic sites with `bcftools v1.13`.

```bash
# Filter out mono-allelic sites.
for CHR in {1..22}; do
bcftools view -m2 -M2 -v snps -Oz -o tgp_altai_merged_filtered_biallelic_chr${CHR}.vcf.gz tgp_altai_merged_filtered_all_sites_chr${CHR}.vcf.gz
done
```

## Trio Results

Run `calc_human_trios.py` to generate site pattern results for every non-African individual.

```bash
# Calculate site pattern counts for all trios.
for CHR in {1..22}; do
python calc_human_trios.py ${CHR}
done
```

## Population Results

To save space I first convert the VCF files to zarr arrays.

```bash
# Convert the VCF files to zarr arrays.
for CHR in {1..22}; do
python vcf_to_zarr.py ${CHR} > ./zarr_arrays/chr${CHR}_altai_zarr_report.txt
done
```

Run `human_site_patterns.py` to generate site pattern results for every non-African population.

```bash
# Calculate site pattern counts for all populations.
python human_site_patterns.py
```

Run `human_bootstrap_replicate.py` to generate a single bootstrapped replicate.

```bash
# Perform bootstrapping.
for REP in {1..1000}; do
python human_bootstrap_replicate.py ${REP}
done
```

## Analysis

A walkthrough of my analysis can be viewed in the `human_analyses_v02.ipynb` notebook.

