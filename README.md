ecre_cds_analysis
=================
Chip-based Library of E. coli Regulatory and Coding Elements 

### Testing

I am testing the github commit feature. Is my key right?

### Prelude: Load RNASeq & DNASeq Data

In order to cut down on the number of bad reads that we are aligning and 
processing, we use what we know about the library sequences to throw away
sequences that are too short. 

First we use `calc_rna_dna.R` to get RNA and DNA read counts. We use a shell
script that wraps `Bowtie` in order to grab the reads:

We throw away any reads that are shorter than the shortest RBS + CDS sequence,
and then calculations to get the RNA/DNA ratio, counts, etc.
