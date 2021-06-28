#######################
#### ORIGINAL CODE ####
#######################

# Original code obtained from https://github.com/DiuTTNguyen/MSI2_HyperTRIBE_codes/tree/master/Pipeline%20code_Yuhueng (tribe_alignment_human.sh)

# The bsub command is used to submit a batch script (see meaning of options here: https://www.ibm.com/docs/en/spectrum-lsf/10.1.0?topic=bsub-options)

#BSUB-J MOLM13_alignment               # name of the job
#BSUB-W 3:00                   	# Time limit in minutes
#BSUB-n 4
#BSUB-R "span[ptile=4]"
#BSUB-M 10

#BSUB-o stdout
#BSUB-o stderr

# Step 1: STAR alignment
ls -d Sample*[0-9] | xargs -I {} sh -c "STAR --genomeLoad NoSharedMemory --genomeDir hg19_star/ --readFilesIn {}/*R1_001.fastq.gz {}/*R2_001.fastq.gz --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outFileNamePrefix {}/ --outStd Log --readFilesCommand zcat"

# ls -d: It is used to display only subdirectories.
# Bracketed characters ([ ]) – matches any occurrence of character enclosed in the square brackets.
# The xargs command is used in a UNIX shell to convert input from standard input into arguments to a command
# The -I option allows you to get the output into a placeholder, and then you can do various things.
# sh is a command language interpreter that executes commands read from a command line string, the standard input, or a specified file. If the -c option is present, then commands are read from string. 

# Meaning of commands :
# --genomeLoad NoSharedMemory : no shared memory usage for genome files
# --genomeDir hg19_star/ : specifies the path to the directory where genome indices are stored
# --readFilesIn {}/*R1_001.fastq.gz {}/*R2_001.fastq.gz : names of the files containing sequences to be mapped (fastq files). If using PE reads, read1 and read2 files have to be supplied
# --runThreadN 4 : defines number of threads to be used for genome generation. Has to be set to the number of available cores on server. 
# --alignIntronMin 70 : minimum intron length (default 20)
# --alignIntronMax 100000 : maximum intron length (default 1000000)
# --outSAMtype BAM SortedByCoordinate : outputs an alignment file that is sorted by coordinates (similar to samtools sort)
# --outFilterMultimapNmax 1 : maximum number of loci the read is allowed to map to. Alignments will be output only if the reads maps to no more loci than this value (default 10)
# --outFilterMultimapScoreRange 0 : the score range below the maximum score for multimapping alignments (default 1)
# --outFilterMismatchNmax 5 : alignment will be output only if it has no more mismatches than this value (default 10)
# --outFileNamePrefix {}/ : output files name prefix
# --outStd Log : determines which output will be directed to standard output
# --readFilesCommand zcat : used to read in compressed (i.e. gz files)

# Step 2: Rename files
ls -d Sample*[0-9] | xargs -I {} sh -c "mv {}/*bam {}.bam"

# Step 3: Index bam files
# The -n option lets you tell xargs to perform one iteration at a time
# samtools index : Index a coordinate-sorted BGZIP-compressed SAM, BAM or CRAM file for fast random access
ls Sample*[0-9].bam | xargs -n 1 samtools index

#################
#### TEST RUN ###
#################

## STAR manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

# Download the data: fasta genome sequence and gtf annotation file
## download hg19 genome: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/ (hg19.fa.gz)
## download hg19 annotation: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/  (hg19.ensGene.gtf.gz)
## unzip the .gz files using gunzip <file>

## 1. Indexing hg19 genome ~1hr
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir hg19-index/ --genomeFastaFiles hg19/hg19.fa --sjdbGTFfile hg19/hg19.ensGene.gtf --sjdbOverhang 49

# --runMode genomeGenerate : directs STAR to run genome indices generation job
# --sjdbGTFfile hg19/hg19.ensGene.gtf : (splice junction database) path to GTF file with annotations
# --sjdbOverhang 49 : length of donor/acceptor sequence on each side of the junctions. Ideally: (mate length - 1). For instance, for Illumina 2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the ideal value is max(ReadLength)-1

## 2. Aligning

# output files will be named like this: Molm13_MIG_R1_Aligned.sortedByCoord.out.bam 
# BAM format explained: https://support.illumina.com/help/BS_App_RNASeq_Alignment_OLH_1000000006112/Content/Source/Informatics/BAM-Format.htm

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_MIG_R1_1.fastq fastq-files/Molm13_MIG_R1_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_MIG_R1_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_MIG_R2_1.fastq fastq-files/Molm13_MIG_R2_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_MIG_R2_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_MIG_R3_1.fastq fastq-files/Molm13_MIG_R3_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_MIG_R3_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_Hyper_dADAR_DCD_R1_1.fastq fastq-files/Molm13_Hyper_dADAR_DCD_R1_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_Hyper_dADAR_DCD_R1_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_Hyper_dADAR_DCD_R2_1.fastq fastq-files/Molm13_Hyper_dADAR_DCD_R2_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_Hyper_dADAR_DCD_R2_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_Hyper_dADAR_DCD_R3_1.fastq fastq-files/Molm13_Hyper_dADAR_DCD_R3_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_Hyper_dADAR_DCD_R3_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_Hyper_dADAR_R1_1.fastq fastq-files/Molm13_Hyper_dADAR_R1_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_Hyper_dADAR_R1_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_Hyper_dADAR_R2_1.fastq fastq-files/Molm13_Hyper_dADAR_R2_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_Hyper_dADAR_R2_

STAR --genomeLoad NoSharedMemory --genomeDir hg19-index/ --readFilesIn fastq-files/Molm13_Hyper_dADAR_R3_1.fastq fastq-files/Molm13_Hyper_dADAR_R3_2.fastq --runThreadN 4 --alignIntronMin 70 --alignIntronMax 100000 --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 1 --outFilterMultimapScoreRange 0 --outFilterMismatchNmax 5 --outStd Log --outFileNamePrefix mapped-reads-OUT/Molm13_Hyper_dADAR_R3_
