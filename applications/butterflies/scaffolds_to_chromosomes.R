#!/usr/bin/Rscript

# Function from Martin et al 2013 to convert scaffold positions to chromosome positions.
as.chromosomes <- function(table,agp) {
    # Initializing empty lists the length of the table.
    chromosome <- vector(length = length(table[,1]))
    chrom_pos <- vector(length = length(table[,1]))
    # For every position in the frequency table.
    for (x in 1:length(table[,1])){
        # Grab the scaffold for the frequency table.
        scaf <- table$scaffold[x]
        # Check if the scaffold is in the mapping file.
        if (table$scaffold[x] %in% agp$scaffold){
            # Identify where the scaffold index in the mapping file.
            y <- which(agp$scaffold == scaf)
            # Set the contig to its corresponding chromosome.
            chromosome[x] <- agp$chromosome[y]
            # If the orientation is positive.
            if (agp$ori[y] == '+'){
                # Add the the scaffold starting position to the contig position.
                chrom_pos[x] <- table$position[x] + agp$start[y]
            }
            # If the orientation is reversed.
            else {
                # Subtract the contig position from the scaffold end position.
                chrom_pos[x] <- agp$end[y] - table$position[x]
            }
        }
        # If not in the mapping file set to NA.
        else{
            chromosome[x] <- NA
            chrom_pos[x] <- NA
        }
        # Update progress.
        if (x %% 100000 == 0){
            print(x)
        }
    }
    # Concatenate the new columns to the frequency table.
    WG_table <- cbind(chromosome,chrom_pos,table)
}


# Read agp for chromosomal info.
agp <- read.delim('./Hmel1-1_hox_RAD_matepair_chromosomes_Zupdated.agp', as.is = T, sep = '\t', header = F)
names(agp) = c('chromosome', 'start', 'end', 'order', 'DN', 'scaffold', 'one', 'length', 'ori')


# Read derived frequency information.
freq_table <- read.csv('./set31.Zupdated.autoANDchrZ.ALLSHARED.SNP.derFreqs.csv.gz')


# Rearrange lines into chromosomes, discarding unmapped scaffolds.
WG_table <- as.chromosomes(freq_table,agp)


# Define the chromosomes of interest.
chromNames <-  c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chrZ')


# Filter out the unmapped chromosomes.
WG_table_subset <- WG_table[WG_table$chromosome %in% chromNames,]


# Factor chromosomes.
WG_table_subset$chromosome <- factor(WG_table_subset$chromosome, levels = chromNames)


# Sort by chromosome first then by position.
WG_table_subset_sorted <- WG_table_subset[
  with(WG_table_subset, order(WG_table_subset$chromosome, WG_table_subset$chrom_pos)),
]


# Filter out sites that don't have derived frequency information.
WG_table_subset_sorted_filtered <- na.omit(WG_table_subset_sorted)


# Output the converted frequency table to a csv.
write.csv(WG_table_subset_sorted_filtered, file = './butterflies_wgs_filtered_der_freqs.csv', row.names = F, quote = F)