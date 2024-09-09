#!/bin/bash

#######################
#### CONFIGURATION ####
#######################

# Number of threads
THREADS=16

# Memory allocation (in GB)
MEMORY=256

# Paths to genome index and input/output directories
GENOME_DIR="/home/linxy29/data/reference/STAR_ref_GRCh38_and_H1N1"
OUTPUT_DIR="extract_RNAedit/"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Function to check if a command succeeded
check_success() {
    if [ $? -ne 0 ]; then
        echo "Error: $1 failed. Exiting."
        exit 1
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

# Get a list of unique sample names by stripping _R1/_R2 from the fastq file names
samples=($(ls *.fastq.gz | sed -E 's/(_R[12]_001\.fastq\.gz)$//' | sort | uniq))

# Number of samples to process
num_samples=${#samples[@]}
echo "Number of samples to process: $num_samples"

#######################
#### ALIGNMENT LOOP ###
#######################

# Loop through each sample and process R1 and R2
for sample in "${samples[@]}"
do
    star_outfile="${OUTPUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam"

    # Check if this sample has already been processed
    if file_exists $star_outfile "STAR alignment"; then
        continue
    fi

    # Define paths for R1 and R2 fastq files
    fastq_R1="${sample}_R1_001.fastq.gz"
    fastq_R2="${sample}_R2_001.fastq.gz"
    
    # Check if both R1 and R2 exist
    if [[ -f $fastq_R1 && -f $fastq_R2 ]]; then
        echo "Processing sample: $sample"
        echo "R1: $fastq_R1"
        echo "R2: $fastq_R2"

        # Run STAR alignment for the sample
        STAR --genomeDir $GENOME_DIR \
             --readFilesIn $fastq_R1 $fastq_R2 \
             --runThreadN $THREADS \
             --outSAMtype BAM SortedByCoordinate \
             --outFilterMultimapNmax 1 \
             --outFileNamePrefix ${OUTPUT_DIR}/${sample}_ \
             --readFilesCommand zcat
        
        # Check for success of the STAR command
        check_success "STAR alignment for $sample"

    else
        echo "Missing R1 or R2 file for sample: $sample"
    fi
done
