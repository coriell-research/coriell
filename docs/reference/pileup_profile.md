# Return the pileup of reads over a set of GRanges

This function computes pileups for BAM files over a set of GRanges. It
returns pileup stats relative to the 5' to 3' direction taking
strandedness of the input range into account.

## Usage

``` r
pileup_profile(
  bamfile,
  x,
  scan_bam_flag,
  max_depth = 10000,
  min_base_quality = 13,
  min_mapq = 1,
  min_nucleotide_depth = 0,
  summarize = TRUE
)
```

## Arguments

- bamfile:

  Path to an indexed BAM file

- x:

  GRanges to calculate pileups over

- scan_bam_flag:

  An instance of Rsamtools::scanBamFlag() used to to create a parameter
  object influencing what fields and which records are imported from a
  (binary) BAM file

- max_depth:

  maximum number of overlapping alignments considered for each position
  in the pileup. See Rsamtools::PileupParam()

- min_base_quality:

  minimum ‘QUAL’ value for each nucleotide in an alignment. Use
  phred2ASCIIOffset to help translate numeric or character values to
  these offsets.

- min_mapq:

  minimum ‘MAPQ’ value for an alignment to be included in pileup.

- min_nucleotide_depth:

  minimum count of each nucleotide (independent of other nucleotides) at
  a given position required for said nucleotide to appear in the result.

- summarize:

  Should the summary stats of the pileups at every relative position
  over all ranges be returned? Default TRUE. If FALSE then a data.table
  containing all measured ranges and relative positions is returned.

## Value

data.table

## Details

For pileup computation, ranges are split by positive and negative strand
and positions in the final data.table are returned relative to the 5'
start in the input ranges. If unstranded ranges are supplied then they
are treated as positively stranded ranges.

## Examples

``` r
if (FALSE) { # \dontrun{

# Select GRanges to compute pileups over
gtf <- rtracklayer::import("/path/to/annotation.gtf")
genes <- gtf[gtf$type == "gene" & gtf$gene_biotype == "protein_coding", ]

# Define parameters for reading in BAM files
flags <- Rsamtools::scanBamFlag(
  isPaired = TRUE,
  isProperPair = TRUE,
  isUnmappedQuery = FALSE,
  hasUnmappedMate = FALSE
  )

# Compute pileups over all gene ranges  
result <- pileup_profile("/path/to/sorted.bam", x = genes, scan_bam_flag = flags)

} # }
```
