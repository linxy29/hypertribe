library(VariantAnnotation)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BiocParallel)
library(parallel)
library(GenomicAlignments)
library(Rsamtools)
library(reshape2)
library(openxlsx)

bam.files <- Sys.glob("*dedupped.bam") # Function to do wildcard expansion (also known as globbing) on file paths. Returns a character vector of matched file paths
print(bam.files)

## Step 1: Filter the SNPs we are interested in

# Input: Given a VCF file from previous step, a list of exons (GRanges object with ranges of all exons within UCSC hg38), and hg38 genome
# Output: This function returns a GRanges object, where each entry is a A/G or T/C SNP mutation that occurs within an exon.
filter.vcf <- function( vcf.file, exons, genome ){
    vcf <- readVcf( vcf.file, genome ) # Read in Variant Call Format (VCF) files. vcf is a dataframe
# Filter all SNPs found in dbSNP
    vcf <- vcf[ !info(vcf)$DB ] #$DB is a logical field (TRUE or FALSE). ! changes the FALSE value to TRUE. This grabs all the novel SNPs, i.e. removes those found in dbSNP
# Only keep A/G mutations in forward strand and T/C mutations in reverse strand
    mask1 <- rowRanges( vcf )$REF == 'A' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'G' %in% seq ) ) # sapply() function in R Language takes list, vector or data frame as input and gives output in vector or matrix. It is useful for operations on list objects and returns a list object of same length of original set. Syntax: sapply(X, FUN
    mask2 <- rowRanges( vcf )$REF == 'T' & sapply( rowRanges( vcf )$ALT, function( seq ) return( 'C' %in% seq ) )
    gr1 <- rowRanges( vcf )[ mask1 ]
    gr2 <- rowRanges( vcf )[ mask2 ]
    strand(gr1) <- '+' # why hard code strand?
    strand(gr2) <- '-'
    a2g.snp <- c( gr1, gr2 ) # combine the 2 Granges objects. I believe these 2 objects are equivalent? Since an A/G transition will lead to a T/C transition on the reverse strand
# Only keep A/G mutations within exons
    a2g.snp <- a2g.snp[ countOverlaps( a2g.snp, exons ) > 0 ]
    return(a2g.snp)
}

vcf.files <- Sys.glob( "*FinalR.vcf" )
# Load all the exons from UCSC gene models
exbygene <- exonsBy( TxDb.Hsapiens.UCSC.hg38.knownGene, "gene" ) # Returns a large compressed GRanges list of all genes within hg38. Grouped by genes.
l <- mclapply( vcf.files, filter.vcf, exbygene, "hg38", mc.cores = 1 ) #exbygene and hg38 are arguments into filter.vcf function
# mclapply is a parallelized version of lapply, it returns a list of the same length as X, each element of which is the result of applying FUN to the corresponding element of X. It relies on forking and hence is not available on Windows unless mc.cores = 1. 
# quite long ~30min

# all A to G and T to C mutations across all vcfs
names(l) <- basename( vcf.files )

# make into a list of GRanges
grl <- GRangesList(l) 
grl
# Rds files store a single R object. According to R documentation: These functions provide the means to save a single R object to a connection (typically a file) and to restore the object, quite possibly under a different name.
# Clicking onto this file automatically generates this R object, so you don't have to rerun all the functions required to generate this object
saveRDS( grl, file = "a2g_snp_filtered.rds" )

# I want unique SNPs only, since there is a possibility of repeated SNPs across the 9 vcfs. Sorted by coordinate
a2g.snp <- sort( unique( unlist( grl, use.names = F ) ) ) 


## Step 2: Counting filtered SNPs across the 9 datasets

# Take bam.files, which is a list of all bam file names, and input into function(bf). The bam files used here are dedupped.bam.
# pileup.res output: "piling up" all the reads at each SNP (across all the bam files) we have obtained previously, into a count.
pileup.res <- mclapply( bam.files, function( bf ){
    bf <- BamFile( bf ) # Use BamFile() to create a reference to a BAM file (and optionally its index). The reference remains open across calls to methods, avoiding costly index re-loading.
    # pileup visits each position in the BAM file, subject to constraints implied by which and flag of scanBamParam. For a given position, all reads overlapping the position that are consistent with constraints dictated by flag of scanBamParam and pileupParam are used for counting. Having "which" indicates to only count the reads for those positions in a2g.snp. Additional flags: we want reads that have mapped mates, is part of a proper pair and are not duplicates. We don't care which strand the read is counted from. Max_depth = maximum number of overlapping alignments considered for each position in the pileup. Min_Base_quality = minimum ‘QUAL’ value for each nucleotide in an alignment. 
    res <- pileup( bf, a2g.snp, 
        scanBamParam=ScanBamParam( flag = scanBamFlag( hasUnmappedMate=F, isProperPair=T, isDuplicate=F ), which = a2g.snp ), 
        pileupParam=PileupParam( distinguish_strands=F, min_base_quality=10, max_depth=1e4 ) )
    return( res ) }, mc.cores = 1 )

saveRDS( pileup.res, file = "pileup.res_bf.rds" )

class(pileup.res) # a list
#pileup.res[[9]] # ninth bam file

names( pileup.res ) <- basename( bam.files ) # renames according to file name

# Saving pileup object for easier future access
#saveRDS( pileup.res, file = "pileup_res_DCD_only.rds" )
saveRDS( pileup.res, file = "pileup_res.rds" ) # use this one

# if( F ){
pileup.res <- lapply( pileup.res, function( df ){
    df <- dcast( df, which_label ~ nucleotide, value.var = 'count', fill = 0, drop = F ) #recasting the df with x variable: which_label and y variable: nucleotide (we want to see how many counts for each which_label). value.var: Name of the column whose values will be filled to cast. fill: Value with which to fill missing cells. drop: FALSE will cast by including all missing combinations.
    df$strand <- as.character( strand( a2g.snp ) ) # converting strand into character and inputting into a column. C/T is on reverse strand, A/G is on forward (positive) strand
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) ) # s: string, d: integer (https://www.gastonsanchez.com/r4strings/c-style-formatting.html)
    stopifnot( all( snp.id == df$which_label ) ) # is TRUE. Means that all snp.ids are accounted for in the which_label of df
    rownames( df ) <- df$which_label # label the SNPs rownames based on their position
    df <- split( df, df$strand ) # split the SNPs based on strand (either forward or reverse)
    df$`+` <- data.frame( ref = 'A', alt = 'G', ref.count = df$`+`$A, alt.count = df$`+`$G, row.names = rownames( df$`+` ) ) # manipulating the positive strand df to obtain counts of each SNP
    df$`-` <- data.frame( ref = 'T', alt = 'C', ref.count = df$`-`$T, alt.count = df$`-`$C, row.names = rownames( df$`-` ) ) # manipulating reverse strand df
    df <- rbind( df$`+`, df$`-` )[ snp.id, ] # combines the positive and negative strand df into one big df, based on their snp.id
    return( df )
} )

head(pileup.res) # each sample is appended with SNP id, ref allele, alt allele, ref.count, alt.count

### Step 3: Annotating SNPs

# Most simplistic allele annotation - just the gene where the SNP is from 
idx <- which( countOverlaps( a2g.snp, exbygene ) == 1 ) # find out which SNP has an overlap with gene annotation. TRUE if exists (overlaps), FALSE if not. The which() function in R returns the position or the index of the value which satisfies the given condition.
a2g.snp.subset <- a2g.snp[idx] # only getting the SNPs that are annotated with gene
ol  <- findOverlaps( a2g.snp.subset, exbygene ) # Finds range overlaps between a GAlignments, GAlignmentPairs, or GAlignmentsList object, and another range-based object. When the query and the subject are GRanges or GRangesList objects, findOverlaps uses the triplet (sequence name, range, strand) to determine which features (see paragraph below for the definition of feature) from the query overlap which features in the subject, where a strand value of "*" is treated as occurring on both the "+" and "-" strand. An overlap is recorded when a feature in the query and a feature in the subject have the same sequence name, have a compatible pairing of strands (e.g. "+"/"+", "-"/"-", "*"/"+", "*"/"-", etc.), and satisfy the interval overlap requirements. 
# Hits object: The as.data.frame method coerces a Hits object to a two column data.frame with one row for each hit, where the value in the first column is the index of an element in the query and the value in the second column is the index of an element in the subject.
entrez.id <- names( exbygene )[ subjectHits(ol) ] # Entrez Gene provides unique integer identifiers for genes and other loci (such as officially named mapped markers) for a subset of model organisms. (Entrez ID is maintained by NCBI)
gene.symbol <- mget( entrez.id, org.Hs.egSYMBOL, ifnotfound = NA ) # get the gene symbol for entrez.id
stopifnot( all( elementNROWS(gene.symbol) == 1 ) ) # each entrez id maps back to only 1 unique gene symbol
gene.symbol <- unlist( gene.symbol ) # becomes a character
anno <- DataFrame( entrez.id, gene.symbol ) # create a dataframe with column1 = entrez.id and column2 = gene.symbol

anno$annotation <- "cds" # add a column called "annotation" and fill the whole column with "cds"
utr3bytx <- threeUTRsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene) #  Get the 3' UTRs grouped by transcript
utr5bytx <- fiveUTRsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene) #  Get the 5' UTRs grouped by transcript
anno$annotation[countOverlaps( a2g.snp.subset, utr5bytx ) > 0] <- 'utr5' # if the snp overlaps with UTR region, change the annotation column to UTR instead of CDS
anno$annotation[countOverlaps( a2g.snp.subset, utr3bytx ) > 0] <- 'utr3'

# Should we go back and check the edits found in dbSNP? Maybe some of them are real after all?
# Add back dbSNP annotation
# dbsnp.gr <- endoapply( grl, function( gr ) gr[ grep( '^rs', names(gr) )] )
# dbsnp.gr <- sort( unique( unlist( dbsnp.gr, use.names = F ) ) )

# anno$dbsnp <- ""
# ol <- findOverlaps( a2g.snp.subset, dbsnp.gr )
# anno$dbsnp[ queryHits( ol ) ] <- names( dbsnp.gr )[ subjectHits( ol ) ]

# a2g.snp.subset: only getting the SNPs that are annotated with gene
mcols( a2g.snp.subset ) <- anno # for each SNP, adding 3 metadata columns from anno (entrez.id, gene.symbol & annotation)


pileup.res.test <- lapply( pileup.res, function( df, a2g.snp ){#l in lapply() stands for list. The difference between lapply() and apply() lies between the output return. The output of lapply() is a list. lapply() can be used for other objects like data frames and lists.
    snp.id <- sprintf( "%s:%d-%d", seqnames( a2g.snp ), start( a2g.snp ), start( a2g.snp ) ) # eg snp.id: chr1:2116555-2116555
    df <- df[ snp.id, ]
    df <- cbind( as.data.frame( a2g.snp), df )
    df$end <- NULL
    df$width <- NULL
    return( df )
}, a2g.snp.subset )
head(pileup.res$Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam)
head(pileup.res.test$Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam) # matched based on snp.id. each SNP is annotated with seqname, start, strand, entrez.id, gene.symbol, annotation (cds, utr)

anno.df <- pileup.res.test[[1]][(1:(ncol(pileup.res.test[[1]])-2))]
pileup.res.test2 <- lapply( pileup.res.test, function( df ) df[-(1:(ncol(pileup.res.test[[1]])-2))] )
pileup.res.test2 <- do.call( "cbind", pileup.res.test2 )
pileup.res.test2 <- cbind( anno.df, pileup.res.test2 )

write.xlsx( pileup.res.test2, file = "snp_counts_dedupped.xlsx" )

# }
