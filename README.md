# hypertribe
Hypertribe run on Molm13 datasets

Original code: https://github.com/DiuTTNguyen/MSI2_HyperTRIBE_codes

## Workflow
1. Align reads to reference genome
+ Tool: STAR (splice-aware aligner)
+ Output: BAM alignment files
2. Prepare BAM files for variant calling
+ Tool: Picard (manipulate HTS data)
+ Add read groups
+ Mark duplicates
3. Perform variant calling on prepped BAM files (all SNPs)
+ Tool: GATK (Genome Analysis ToolKit for variant discovery)
+ Output: VCF (variant caller format) files
4. Filter SNPs and read counting
+ A/G, T/C
