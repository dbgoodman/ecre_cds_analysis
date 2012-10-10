<link href="http://kevinburke.bitbucket.org/markdowncss/markdown.css" rel="stylesheet"></link>



# Step 02

## Motivation

In [Step 02](02_load_dna_rna.html), we aligned all the DNA and RNA reads directly to a reference library. One of the take home messages from that effort was that there are a lot of mismatching reads, and a lot of them are likely due to synthesis, and not sequencing. 

We are going to compute properties of each construct based on its nucleotide and DNA sequence, and so we are more intrested here in library diversity and not the identities of the sequences themselves. For instance, we don't care so much that the CDS sequence consists precisely the same amino acids as the front of the BamA gene, so much as that it has X secondary structure, Y codon adaptation index, etc. 

I'm going to (quickly) check to see if it is feasible to see if we can proceed without aligning both the RNA and DNA directly to the reference library. Perhaps instead we can first pull out all DNA with sufficient reads, and then subsequently align the RNA to those new 'reference' reads.

### How do we choose a DNA read count cutoff, if we don't align to references?

In order to observe low RNA:DNA ratios we need to limit ourselves to DNA constructs with many DNA reads. To be sure we are looking at a synthesis and not a sequencing error in the RNA, we should see at least 5 reads in both replicates. 

Looking at the `plot_rna_dna.rna_ratio.tech_replicate.pdf` from the 202 analysis, it looks like a good ratio cutoff where most reads fall is about 2 `log10(-1.3)`, or if we have 5 read, then approximately 100 DNA reads per construct. If we don't see at least 5, then we can say the ratio is below 1:20. 

This is pretty much what we did on the 202 analysis, except in that case, almost all of the reads with low DNA occurred because they were high expressing constructs, causing impaired cell growth. This meant that for most of those sequences, the RNA would also be high, so we would not run into this problem often. 

### How many reference constructs have at least 100 reads per replicate?

We can of course see how many constructs that are in the `dna.subset` have over 100 reads:



```r
# Number of contigs matching perfectly with >=100 reads in both
# replicates:
sum(with(dna.subset, Mismatches.len == 0 & Count.A >= 100 & Count.B >= 
    100))
```



```
## [1] 13922
```



```r

# Number of constructs in library:
dim(lib_seqs)[1]
```



```
## [1] 14234
```



```r

# Designed sequences lost by this approach:
dim(lib_seqs)[1] - sum(with(dna.subset, Mismatches.len == 0 & Count.A >= 
    100 & Count.B >= 100)) - length(dna.missing)
```



```
## [1] 297
```




So we'd lose almost 300 of the 14,234 sequences due to low DNA counts. 

### How many new constructs might we gain?

How many would we gain? Let's take a look at the raw counts file. Under this scheme, we want to discard any reads with Ns, since we won't be able to assign them easily. We could default them to the designed reference construct, but that would be a pain. We also are setting size restrictions between 70 and 100 bases. 

```bash
lib_prefix=/scratch/dbg/ecre/fa/203.norestrict
out_prefix=/scratch/dbg/ecre/203_hs
dna_prefix=/scratch/dbg/ecre/ct/203_hsdna/203_hsdna.counts
rna_prefix=/scratch/dbg/ecre/ct/203_hsrna/203_hsrna.counts

#Count unique contigs with at least 100 reads in both replicates
less -S $dna_prefix.fa | perl -pe 'chomp; s/>/\n>/; 
        s/^([ATGC])/\t$1/;' \
    | perl -ne '@l = split; !(/N/) && $l[2] >= 100 && $l[3] >= 100 
        && length($l[4]) < 100 && length($l[4]) > 70 && print $_;' \
    | sort -nrk2 | cut -f2 | uniq -c \
    | cut -f1 | perl -nle '$sum += $_ } END { print $sum'
    
#Count unique contigs with at least 100 reads in both replicates
less -S $dna_prefix.fa | perl -pe 'chomp; s/>/\n>/; 
        s/^([ATGC])/\t$1/;' \
    | perl -ne '@l = split; !(/N/) && $l[2] >= 100 && $l[3] >= 100 
        && length($l[4]) < 100 && length($l[4]) > 70 && print $_;' \
    | sort -nrk2 | cut -f2 | uniq -c \
    | cut -f1 | perl -nle '$sum += $_ } END { print $sum'

#Output: 16972

#Count unique contigs with at least 50 reads in both replicates
less -S $dna_prefix.fa | perl -pe 'chomp; s/>/\n>/; 
        s/^([ATGC])/\t$1/;' \
    | perl -ne '@l = split; !(/N/) && $l[2] >= 50 && $l[3] >= 50 
        && length($l[4]) < 100 && length($l[4]) > 70 && print $_;' \
    | sort -nrk2 | cut -f2 | uniq -c \
    | cut -f1 | perl -nle '$sum += $_ } END { print $sum'

#Output: 27154
```

So we end up with 16972 sequences with 100+ DNA reads. That's an extra 2,738 reads. If we lower our cutoff to 50 (which are still likely synthesis errors, but we'd have a smaller RNA ratio dynamic range), then that number goes up to 27154, or double our library size. 

### Caveats to this approach

Now, there are some caveats here. We can't use all of these sequences for a few reasons:

1. We should probably limit ourselves to sequences with perfect promoters. I don't think this won't hurt us too much, because the mismatch distribution histograms from Step 02 showed that not too many of the mismatches occur before base 40. 
2. The CDS need to be in frame, from the first ATG onward. This is easy enough to check after splitting at the ATG. I'm not sure whether I should throw away or keep shorter CDSs; it will probably depend on how many there are. 
3. Dealing with multiple potential start sites might be problematic. For instance, if a mutation creates an extra ATG somewhere else close to the Shine-Dalrarno, it would complicate the amino acid analysis. To ensure that the start site does not change, I should remove sequences that create new ATGs close to the RBS in different frames, and make sure that none of the codons in the CDS were mutated to ATG. 
4. I also want to throw away any sequences with stop codons in them. 
5. Some of these sequences (42 + mismatches and contaminants) are the 202-spike-in set.
6. As Sri and I discussed earlier, we are hiding any errors that are not part of the variable region, like the fluorescent coding sequences on the backbone.

>I'll show this to Sri. I think we could probably do it this way. It resolves all of the problems that we've seen with the DNA mismatch sequences and the extra sequences would likely enhance the power of our analyses.
