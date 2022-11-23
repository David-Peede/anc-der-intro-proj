# Canid Application

This directory contains all the code to replicate the canid analyses from Peede et al. 202X, which I will now outline.

## Download

The canid data is publicly available and can be downloaded from the [dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.sk3p7) associated with [Freedman et al. 2014](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004016).

## Filter

I used `bcftools v1.13` to filter out sites that were not bi-allelic and sites that did not pass the genomic filter specified in [Freedman et al. 2014](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004016).

```bash
# Filter out sites that do not pass the genomic feature filter.
for CHR in {1..38}; do
bcftools view -M2 -i 'GF = 1' -Oz -o canid_merged_filtered_chr${CHR}.vcf.gz simp_2012.5.18_6canid_merged_chr${CHR}_cf31.vcf.gz
done
```

## Site Patterns

I used `canid_site_patterns.py` to compute the observed results.

```bash
# Calculate site patterns.
python canid_site_patterns.py
```

## Bootstrapping

To help speed up the bootstrapping process I first create a `.csv` per chromosome with all the positions that based quality control.

```bash
# Create csvs with position info per chromosome.
python canid_vcf_positions.py
```

Run `canid_bootstrap_replicate.py` to generate a single bootstrapped replicate.

```bash
# Perform bootstrapping
for REP in {1..1000}; do
python canid_bootstrap_replicate.py ${REP}
done
```

## Analysis

A walkthrough of my analysis can be viewed in the `canid_analyses.ipynb` notebook.

