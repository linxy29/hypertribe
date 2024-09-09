#!/bin/bash

# Paths
REFERENCE="/home/linxy29/data/reference/GRCh38_and_H1N1.fa"
DBSNP="/home/linxy29/data/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
THREADS=30

# Function to check if a command succeeded
check_success() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed. Exiting."
        exit 1
    else
        echo "$1 completed successfully."
    fi
}

# Function to check if a file exists
file_exists() {
    if [ -f $1 ]; then
        echo "Skipping $2, file $1 already exists."
        return 0
    else
        return 1
    fi
}

# Iterate over all *_Aligned.sortedByCoord.out.bam files
for bam in *_Aligned.sortedByCoord.out.bam; do
    sample=$(basename $bam _Aligned.sortedByCoord.out.bam)

    # Step 1: Add or Replace Read Groups
    output_rg="${sample}_rg_added.bam"
    if ! file_exists $output_rg "AddOrReplaceReadGroups for $sample"; then
        picard AddOrReplaceReadGroups \
            I=$bam \
            O=$output_rg \
            RGID=$sample RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=$sample
        check_success "Adding read groups for $sample"
    fi

    # Step 2: Mark Duplicates
    output_dedup="${sample}_dedupped.bam"
    if ! file_exists $output_dedup "MarkDuplicates for $sample"; then
        picard MarkDuplicates \
            I=$output_rg \
            O=$output_dedup \
            CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=${sample}.metrics
        check_success "Marking duplicates for $sample"
    fi

    # Step 3: SplitNCigarReads
    output_split="${sample}_split.bam"
    if ! file_exists "$output_split" "SplitNCigarReads for $sample"; then
        gatk SplitNCigarReads \
            -R /home/linxy29/data/reference/GRCh38_and_H1N1.fa \
            -I $output_dedup \
            -O $output_split 
        check_success "Splitting reads for $sample"
    fi

    # Step 4: Call SNPs and Indels using HaplotypeCaller
    output_vcf="${sample}.vcf"
    if ! file_exists $output_vcf "HaplotypeCaller for $sample"; then
        gatk HaplotypeCaller \
            -R $REFERENCE \
            -I $output_split \
            --dont-use-soft-clipped-bases TRUE\
            -stand-call-conf 10.0 \
            --dbsnp $DBSNP \
            -O $output_vcf \
            --native-pair-hmm-threads $THREADS
        check_success "Calling variants for $sample"
    fi

    # Step 5: Annotate variants using VariantAnnotator
    output_annotated_vcf="${sample}_annotated.vcf"
    if ! file_exists $output_annotated_vcf "VariantAnnotator for $sample"; then
        gatk VariantAnnotator \
            -R $REFERENCE \
            -V $output_vcf \
            -O $output_annotated_vcf \
            -A FisherStrand \
            -A QualByDepth
        check_success "Annotating variants for $sample"
    fi


    # Step 6: Variant Filtration
    output_filtered_vcf="${sample}_filtered.vcf"
    if ! file_exists $output_filtered_vcf "VariantFiltration for $sample"; then
        gatk VariantFiltration \
            -R $REFERENCE \
            -V $output_annotated_vcf \
            -window 35 -cluster 3 \
            --filter-name "FS" --filter-expression "FS > 30.0" \
            --filter-name "QD" --filter-expression "QD < 2.0" \
            -O $output_filtered_vcf
        check_success "Filtering variants for $sample"
    fi

done

echo "Pipeline completed successfully."
