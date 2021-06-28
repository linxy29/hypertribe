#######################
#### ORIGINAL CODE ####
#######################

#BSUB-J MOLM13_variant_call        # name of the job
#BSUB-W 6:00                   	# Time limit 
#BSUB-n 9
#BSUB-R "span[ptile=9]"
#BSUB-M 30

#BSUB-o stdout
#BSUB-o stderr

# Load Java 8
module add java

# Picard: add RG & mark duplicates
# parallel:  allows the user to execute shell scripts or commands in parallel.
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/picard.jar AddOrReplaceReadGroups I={}.bam O={}.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample'

# Note: this step takes ~30G memory for each file; Be sure to reserve enough RAM
ls -d Sample_*[0-9] | xargs -I {} sh -c 'java -jar ~/jartools/picard.jar MarkDuplicates I={}.rg_added.bam O={}.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics'

############ ############ ############ ignore ############ ############ ############
# hg19 genome we used for alignment is slightly different than the one used by picard 
# For now we have to rename and reorder the chrs to make them consistent
# Someone should try to fix this later
ls -d Sample_*[0-9] | parallel 'samtools reheader header.txt {}.dedupped.bam > {}.reheader.bam' #rename
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/picard.jar ReorderSam I={}.reheader.bam O={}.reordered.bam R=hg19_anno/ucsc.hg19.fasta CREATE_INDEX=TRUE' #reorder
############ ############ ############ ignore ############ ############ ############

# Split reads 
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T SplitNCigarReads -R hg19_anno/ucsc.hg19.fasta -I {}.reordered.bam -o {}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'
# Call SNPs & annotate
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19_anno/ucsc.hg19.fasta -I {}.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o {}.HaplotypesR.vcf --dbsnp hg19_anno/dbsnp_138.hg19.vcf'
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T VariantFiltration -R hg19_anno/ucsc.hg19.fasta -V {}.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o {}.FinalR.vcf'


#################
#### TEST RUN ###
#################

# Pipeline follows this "RNAseq short variant discovery (SNPs + Indels)": https://gatk.broadinstitute.org/hc/en-us/articles/360035531192?id=4067 
# Purpose of variant_calling pipeline: Identify short variants (SNPs and Indels) in RNAseq data.

## Check java version 
java -version 
# should be 1.8 and above 

# Download picard https://broadinstitute.github.io/picard/ 
# The tools, which are all listed further below, are invoked as follows:
# java jvm-args -jar picard.jar PicardToolName OPTION1=value1 OPTION2=value2...

# Step 1: Picard add RG (read group)
# Time taken: 5min each
java -jar ~/jartools/picard.jar AddOrReplaceReadGroups I={}.bam O={}.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

# AddOrReplaceReadGroups (Tool that manipulates read data in SAM, BAM or CRAM format): Assigns all the reads in a file to a single new read-group.
# The program enables only the wholesale assignment of all the reads in the INPUT to a single read-group. If your file already has reads assigned to multiple read-groups, the original RG value will be lost.
# What is a read group: In the simple case where a single library preparation derived from a single biological sample was run on a single lane of a flow cell, all the reads from that lane run belong to the same read group. When multiplexing is involved, then each subset of reads originating from a separate library run on that lane will constitute a separate read group.
# These tags, when assigned appropriately, allow us to differentiate not only samples, but also various technical features that are associated with artifacts. With this information in hand, we can mitigate the effects of those artifacts during the duplicate marking and base recalibration steps. The GATK requires several read group fields to be present in input files and will fail with errors if this requirement is not satisfied. 
# Read more about RG here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_Hyper_dADAR_DCD_R1_Aligned.sortedByCoord.out.bam O=Molm13_Hyper_dADAR_DCD_R1_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_Hyper_dADAR_DCD_R2_Aligned.sortedByCoord.out.bam O=Molm13_Hyper_dADAR_DCD_R2_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_Hyper_dADAR_DCD_R3_Aligned.sortedByCoord.out.bam O=Molm13_Hyper_dADAR_DCD_R3_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_Hyper_dADAR_R1_Aligned.sortedByCoord.out.bam O=Molm13_Hyper_dADAR_R1_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_Hyper_dADAR_R2_Aligned.sortedByCoord.out.bam O=Molm13_Hyper_dADAR_R2_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_Hyper_dADAR_R3_Aligned.sortedByCoord.out.bam O=Molm13_Hyper_dADAR_R3_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_MIG_R1_Aligned.sortedByCoord.out.bam O=Molm13_MIG_R1_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_MIG_R2_Aligned.sortedByCoord.out.bam O=Molm13_MIG_R2_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample

java -jar picard.jar AddOrReplaceReadGroups I=../alignment/mapped-reads-OUT/Molm13_MIG_R3_Aligned.sortedByCoord.out.bam O=Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.bam RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample


# Step 2: Mark duplicates.
# Time: ~20min each
# Results: 20% duplication
## Identifies duplicate reads. This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA.  Duplicates can arise during sample preparation e.g. library construction using PCR.
## Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.  These duplication artifacts are referred to as optical duplicates.
## Explaining PCR duplication: https://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/ 
## How MarkDuplicates work: https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-

java -jar picard.jar MarkDuplicates I=Molm13_Hyper_dADAR_DCD_R1_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_Hyper_dADAR_DCD_R1_Aligned.sortedByCoord.out.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_Hyper_dADAR_DCD_R1.metrics

# M= : File to write duplication metrics to
# CREATE_INDEX=true : Whether to create a BAM index when writing a coordinate-sorted BAM file.
# VALIDATION_STRINGENCY=SILENT : Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.

java -jar picard.jar MarkDuplicates I=Molm13_Hyper_dADAR_DCD_R2_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_Hyper_dADAR_DCD_R2_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_Hyper_dADAR_DCD_R2.metrics

java -jar picard.jar MarkDuplicates I=Molm13_Hyper_dADAR_DCD_R3_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_Hyper_dADAR_DCD_R3_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_Hyper_dADAR_DCD_R3.metrics

java -jar picard.jar MarkDuplicates I=Molm13_Hyper_dADAR_R1_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_Hyper_dADAR_R1_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_Hyper_dADAR_R1.metrics

java -jar picard.jar MarkDuplicates I=Molm13_Hyper_dADAR_R2_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_Hyper_dADAR_R2_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_Hyper_dADAR_R2.metrics

java -jar picard.jar MarkDuplicates I=Molm13_Hyper_dADAR_R3_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_Hyper_dADAR_R3_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_Hyper_dADAR_R3.metrics

java -jar picard.jar MarkDuplicates I=Molm13_MIG_R1_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_MIG_R1_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_MIG_R1.metrics

java -jar picard.jar MarkDuplicates I=Molm13_MIG_R2_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_MIG_R2_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_MIG_R2.metrics

java -jar picard.jar MarkDuplicates I=Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.bam O=Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=Molm13_MIG_R3.metrics


## Install GATK: https://gatk.broadinstitute.org/hc/en-us/articles/360036194592 https://console.cloud.google.com/storage/browser/gatk-software/package-archive/gatk/?pli=1
tar -xvf package-archive_gatk_GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar

# Step 3: Split reads 
# Time: ~40min

# SplitNCigarReads: Splits reads that contain Ns in their CIGAR string (For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments, the interpretation of N is not defined.)
# Splits reads that contain Ns in their cigar string (e.g. spanning splicing events in RNAseq data). Identifies all N cigar elements and creates k+1 new reads (where k is the number of N cigar elements). The first read includes the bases that are to the left of the first N element, while the part of the read that is to the right of the N (including the Ns) is hard clipped and so on for the rest of the new reads. Used for post-processing RNA reads aligned against the full reference.

# SplitNCigarReads is developed specially for RNAseq, which splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions --> hard-clip = totally remove those sequences. (https://github.com/broadgsa/gatk/blob/master/doc_archive/methods/Calling_variants_in_RNAseq.md)

# At this step we also add one important tweak: we need to reassign mapping qualities, because STAR assigns good alignments a MAPQ of 255 (which technically means “unknown” and is therefore meaningless to GATK). So we use the GATK’s ReassignOneMappingQuality read filter to reassign all good alignments to the default value of 60. This is not ideal, and we hope that in the future RNAseq mappers will emit meaningful quality scores, but in the meantime this is the best we can do. In practice we do this by adding the ReassignOneMappingQuality read filter to the splitter command.
# Devon: "Don't literally use Integer0to255 as the value for --outSAMmapqUnique. Either use a single value, e.g., --outSAMmapqUnique 50, or don't bother with the argument (the default is 255). It's extremely likely that Integer0to255 caused the "unique alignment" MAPQ to be set to 0."
# Finally, be sure to specify that reads with N cigars should be allowed. This is currently still classified as an "unsafe" option, but this classification will change to reflect the fact that this is now a supported option for RNAseq processing.

# Most GATK tools additionally require that the main FASTA file be accompanied by a dictionary file ending in .dict and an index file ending in .fai, because it allows efficient random access to the reference bases.
# Create fasta dict and fai file: https://gatk.broadinstitute.org/hc/en-us/articles/360035531652?id=1601
samtools faidx ~/hypertribe/alignment/hg19/hg19.fa
# Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a header but no SAMRecords, and the header contains only sequence records. https://gatk.broadinstitute.org/hc/en-us/articles/360037422891-CreateSequenceDictionary-Picard-
java -jar picard.jar CreateSequenceDictionary R=~/hypertribe/alignment/hg19/hg19.fa O=~/hypertribe/alignment/hg19/hg19.dict

# Actual splitting reads: 
java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R1_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R1_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

# MAPQ: -10log(10^-6) = 60 (MAPQ equals −10 log10 Pr{mapping position is wrong})
# -T,--analysis_type <analysis_type>: Name of the tool to run
# -R,--reference_sequence <reference_sequence>: Reference sequence file
# -rf,--read_filter <read_filter>: Filters to apply to reads before analysis
# -RMQF & -RMQT : from 255 to 60
# ALLOW_N_CIGAR_READS : specify that reads with N cigars should be allowed

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R2_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R2_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R3_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R1_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R1_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R2_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R2_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R3_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_MIG_R1_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_MIG_R1_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_MIG_R2_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_MIG_R2_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.dedupped.bam -o ~/hypertribe/variant_calling/Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS


# Step 4: Call SNPs & annotate 
# Time: ~3-4hrs 
# HaplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
# The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. In other words, whenever the program encounters a region showing signs of variation, it discards the existing mapping information and completely reassembles the reads in that region. This allows the HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when they contain different types of variants close to each other.

# How variants are called:
# 1. Find active regions (region showing signs of variation)
# 2. Re-assemble reads at active regions using De-Bruijn graphs (De Bruijn graphs are used for de novo assembly of sequencing reads into a genome) --> edges with large weights are likely to be variant haplotypes; prune edges that are not well-supported (<2 supporting reads)
# 3. Traverse pruned graph and select best haplotypes (calculate likelihood score of each path)
# 4. Identify potential variant sites by aligning reads against each of a list of plausible haplotypes --> reconstruct a CIGAR string for that haplotype --> produces a list of potential variation sites that will be put through the variant modelling step next


# Acquire reference dbsnp vcf file
# download hg19 dbsnp138 vcf from broad institute google cloud (hg19_v0_Homo_sapiens_assembly19.dbsnp138.vcf): https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false
# don't use the .idx file provided at this link. Create our own since we need to adhere to chromosome naming convention

## adding chr to start of vcf file. Chromosome naming convention
awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' hg19_v0_Homo_sapiens_assembly19.dbsnp138.vcf > hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf

# Compressing VCF files with BGZip and indexing it with Tabix is the standard way VCF files are stored
# bgzip vcf file:
bgzip -c hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf > hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz
# index vcf file:
tabix -fp vcf hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz

# bgzip does blockwise (therefore the b in bgzip) compression of the file, which tabix relies on. That enables tabix to quickly retrieve data by only very partially decompressing a (sometimes huge) file, guided by the index. gzip does not do blockwise compression.

# -stand_call_conf 10.0: The minimum phred-scaled confidence threshold at which variants should be called (default 30)
# -T HaplotypeCaller : Call germline SNPs and indels via local re-assembly of haplotypes (How HaplotypeCaller works: https://gatk.broadinstitute.org/hc/en-us/articles/360056969012-HaplotypeCaller)
# -dontUseSoftClippedBases : https://www.biostars.org/p/109333/ --> explains soft and hard clipped bases. "Hard masked bases do not appear in the SEQ string, soft masked bases do. So, if your cigar is: 10H10M10H then the SEQ will only be 10 bases long. If your cigar is 10S10M10S then the SEQ and base-quals will be 30 bases long. Hard masked bases do not appear in the SEQ string, soft masked bases do. So, if your cigar is: 10H10M10H then the SEQ will only be 10 bases long. If your cigar is 10S10M10S then the SEQ and base-quals will be 30 bases long. In the case of soft-masking, even though the SEQ is present, it is not used by variant callers and not displayed when you view your data in a viewer. In either case, masked bases should not be used in calculating coverage. Both of these maskings are different from deletions. Masking simply means the part of the read can not be aligned to the genome (simplified, but a reasonable assumption for most cases, I think). A deletion means that a stretch of genome is not present in the sample and therefore not in the reads. "

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R1_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R1.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R2_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R2.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_DCD_R3.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R1_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R1.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R2_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R2.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_Hyper_dADAR_R3.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz 

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_MIG_R1_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_MIG_R1.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz 

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_MIG_R2_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_MIG_R2.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz 

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I ~/hypertribe/variant_calling/Molm13_MIG_R3_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam -dontUseSoftClippedBases -stand_call_conf 10.0 -o ~/hypertribe/variant_calling/Molm13_MIG_R3.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz 

# Attempt at parallel: looks like this works
ls -d Molm13_MIG_R*[1-2]_Aligned.sortedByCoord.out.rg_added.dedupped.split.bam | parallel 'java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ~/hypertribe/alignment/hg19/hg19.fa -I {} -dontUseSoftClippedBases -stand_call_conf 10.0 -o {}.HaplotypesR.vcf --dbsnp ~/hypertribe/variant_calling/hg19_v0_Homo_sapiens_assembly19.dbsnp138.chr.vcf.gz'

# Step 5: Variant Filtration
# Input: A VCF of variant calls to filter & One or more filtering expressions and corresponding filter names.
# Output: A filtered VCF in which passing variants are annotated as PASS and failing variants are annotated with the name(s) of the filter(s) they failed. 
# This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.
ls -d Sample_*[0-9] | parallel 'java -jar ~/jartools/GenomeAnalysisTK.jar -T VariantFiltration -R hg19_anno/ucsc.hg19.fasta -V {}.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o {}.FinalR.vcf'


# -filter,--filterExpression <filterExpression>: One or more expression used with INFO fields to filter (Composing filtering expressions can range from very simple to extremely complicated depending on what you're trying to do.)
# -filterName,--filterName <filterName>: Names to use for the list of filters (user-defined)
# -cluster,--clusterSize <clusterSize>: The number of SNPs which make up a cluster. Must be at least 2. Works together with the --cluster-window-size argument.
# -window,--clusterWindowSize <clusterWindowSize>: The window size (in bases) in which to evaluate clustered SNPs. Works together with the --cluster-size argument. To disable the clustered SNP filter, set this value to less than 1.

# Filters applied: 
# 1. FS: ##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias"> (p-value: The p-value gives us the probability of observing the set of results we obtained if the null hypothesis were true, i.e. getting those results purely by chance.) [https://gatk.broadinstitute.org/hc/en-us/articles/360035532152-Fisher-s-Exact-Test]
# - In GATK, we use Fisher’s Exact Test to calculate the FisherStrand annotation, which is an indicator of strand bias, a common source of artifactual calls. The test determines whether there is a difference in the number of reads that support the reference allele and alternate allele on each strand (i.e. number of reads in forward and reverse orientation). The value is reported in the FisherStrand annotation, FS in the VCF.
# - In this example, we want to determine if there is a difference in the number of reads that support the reference allele and alternate allele on each strand. Our null hypothesis is that there is no difference in the number of reads that support the reference allele and alternate allele on each strand (there is no strand bias). The lower the p-value, the less likely we are to believe that there is no strand bias.
# - Converting to Phred scaled p-value: In the GATK context we still want to transform our FS annotation value to Phred scale for convenience before writing it out to the output VCF. To get the Phred-scaled p-value, we simply plug in the p-value of 0.1 into the Phred equation.
# - Note if we had a p-value of 1, meaning there is a 100% chance of there being no bias, the Phred score would be 0. So, a Phred-score closer to 0 means there is a lower chance of there being bias. Higher FS values therefore indicate more bias.

# 2. QD: ##INFO=<ID=QD,Number=1,Type=Float,Description="Variant Confidence/Quality by Depth"> [https://gatk.broadinstitute.org/hc/en-us/articles/360036348492-AS-QualByDepth]
# - Variants for which quality (QUAL) divided by depth of coverage (DP)
# - This annotation puts the variant confidence QUAL score into perspective by normalizing for the amount of coverage available. Because each read contributes a little to the QUAL score, variants in regions with deep coverage can have artificially inflated QUAL scores, giving the impression that the call is supported by more evidence than it really is. To compensate for this, we normalize the variant confidence by depth, which gives us a more objective picture of how well supported the call is.
# - The QD is the QUAL score normalized by allele depth (AD) for a variant. For a single sample, the HaplotypeCaller calculates the QD by taking QUAL/AD. 

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_Hyper_dADAR_DCD_R1.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_Hyper_dADAR_DCD_R1.FinalR.vcf

# I want to REMOVE QD < 2 and FS > 30 (e.g. variant with QD = 0.62 will fail the filter)

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_Hyper_dADAR_DCD_R2.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_Hyper_dADAR_DCD_R2.FinalR.vcf

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_Hyper_dADAR_DCD_R3.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_Hyper_dADAR_DCD_R3.FinalR.vcf

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_Hyper_dADAR_R1.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_Hyper_dADAR_R1.FinalR.vcf

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_Hyper_dADAR_R2.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_Hyper_dADAR_R2.FinalR.vcf

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_Hyper_dADAR_R3.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_Hyper_dADAR_R3.FinalR.vcf

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_MIG_R1.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_MIG_R1.FinalR.vcf

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_MIG_R2.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_MIG_R2.FinalR.vcf

java -jar ~/hypertribe/variant_calling/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T VariantFiltration -R ~/hypertribe/alignment/hg19/hg19.fa -V Molm13_MIG_R3.HaplotypesR.vcf -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o Molm13_MIG_R3.FinalR.vcf

################
#####VCF######## 
################
# VCF file explanation: https://www.youtube.com/watch?v=2P4EItXCtFI&ab_channel=LiquidBrainLiquidBrain
# Comparison btwn reference and sample, to find out base pair differences/variations
# Similar to an MSA, where the first row is the reference sequence. VCF will summarise all the variations within that MSA instead of repeating everything
# Types of genetic variations:
# - transitions, transversions, indels
# - structural variations (insertion, inversion, deletion, duplication, cnv)
# Eg. of when VCF can be used:
# We have 2 sets of variant calls (vs. a reference sequence) and we need to decide which is control and tumour sample. We decide whether the variant frequency is non-zero in tumour but zero in control.

# Structure of VCF file
## 1. Metadata: preceded with a ##
## 2. Header line: single #
## 3. Data: no preceding symbol

## Metadata:
# FILTER: filters that have been applied to the data, and each variant can either pass or fail that filter. The description of the filter is in the metadata and represented by an ID
# INFO: explains the abbreviations in the INFO column

## Data:
# REF: is the base at that position in the reference sequence
