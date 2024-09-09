#!/bin/sh

#----- Update the following variables ----
star_indices="/home/linxy29/data/reference/STAR_ref_GRCh38_and_H1N1"
avgquality="25"
threads=16 # Set the number of threads to use for STAR and SAMtools
#------End update variables -----

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

# Use the loop below to run the code for multiple fastq.gz files
filelist=`ls *.fastq.gz`
for file in ${filelist[@]}
do

# The input FASTQ file should have a .fastq.gz prefix, for example s2_mRNA.fastq.gz or HyperTRIBE_rep1.fastq.gz
#file=$1 uncomment this if run for one file
prefix=${file%.fastq.gz*}
trim_input=$file
trim_outfile=$prefix.trim.fastq 

# Step 1: Trimming low-quality bases using Trimmomatic
if ! file_exists $trim_outfile "Trimmomatic"; then
    echo "Running Trimmomatic..."
    trimmomatic SE -phred33 $trim_input $trim_outfile HEADCROP:6 LEADING:25 TRAILING:25 AVGQUAL:$avgquality MINLEN:19
    check_success "Trimmomatic"
fi

# Step 2: Aligning reads with STAR
star_outfile=$prefix"_Aligned.out.sam"
if ! file_exists $star_outfile "STAR alignment"; then
    input=$trim_outfile
    echo "Running STAR alignment..."
    STAR --runThreadN $threads --outFilterMismatchNoverLmax 0.07 --outFileNamePrefix $prefix"_" --outFilterMatchNmin 16 --outFilterMultimapNmax 1 --genomeDir $star_indices --readFilesIn $input
    check_success "STAR alignment"
fi

# Move the STAR output to the desired SAM filename
output=$prefix".sam"
if ! file_exists $output "Move STAR output"; then
    mv $star_outfile $output
fi

# Step 3: Filtering low-quality alignments with SAMtools
filtered_sam=$prefix"_highquality.sam"
if ! file_exists $filtered_sam "SAMtools view (filtering SAM)"; then
    echo "Filtering low-quality alignments..."
    samtools view -@ $threads -Sh -q 10 $output > $filtered_sam
    check_success "SAMtools view (filtering SAM)"
    mv $filtered_sam $output
fi

# Step 4: Converting SAM to BAM
bam_out=$prefix".bam"
if ! file_exists $bam_out "SAM to BAM conversion"; then
    echo "Converting SAM to BAM..."
    samtools view -@ $threads -bhS $output > $bam_out
    check_success "SAMtools view (SAM to BAM conversion)"
    rm $output
fi

# Step 5: Sorting BAM file before using Picard
sort_out=$prefix".sort.bam"
if ! file_exists $sort_out "SAMtools sort (pre-Picard)"; then
    echo "Sorting BAM file..."
    samtools sort -@ $threads $bam_out -o $sort_out
    check_success "SAMtools sort"
    rm $bam_out
fi

# Step 6: Removing duplicates with Picard
dupremove_bam=$prefix"_nodup.bam"
if ! file_exists $dupremove_bam "Picard MarkDuplicates"; then
    echo "Removing duplicates with Picard..."
    picard MarkDuplicates INPUT=$sort_out OUTPUT=$dupremove_bam METRICS_FILE=dup.txt VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true TMP_DIR=tmp ASSUME_SORTED=true
    check_success "Picard MarkDuplicates"
    rm $sort_out
fi

# Step 7: Sorting the BAM file after duplicate removal
if ! file_exists $sort_out "SAMtools sort (after Picard)"; then
    echo "Sorting BAM file after duplicate removal..."
    samtools sort -@ $threads $dupremove_bam -o $sort_out
    check_success "SAMtools sort (after Picard)"
    rm $dupremove_bam
fi

# Step 8: Creating a sorted SAM file and indexing BAM
sorted_sam=$prefix".sort.sam"
if ! file_exists $sorted_sam "SAMtools view (BAM to SAM conversion)"; then
    echo "Creating sorted SAM file and indexing BAM..."
    samtools view -@ $threads -h $sort_out > $sorted_sam
    check_success "SAMtools view (BAM to SAM conversion)"
fi

if ! file_exists $sort_out.bai "SAMtools index"; then
    samtools index $sort_out
    check_success "SAMtools index"
fi

# Step 9: Completion message
echo "Done with STAR mapping and PCR duplicate removal with Picard."
echo "Created SAM file: $prefix.sort.sam"

done