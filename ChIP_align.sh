#!/bin/bash -l
#
# Set the name of the job
#SBATCH --job-name=ChIPalign
#
# Set the maximum memory allowed
#SBATCH --mem=20G
#
# Set the maximum run time
#SBATCH -t 48:00:00
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email
#
# The number of threads we will require
#SBATCH -n 8
#
# Set output and error log files
#SBATCH -o /<path_to_working_directory>/step1/Log/logOutAlignAll_ChIP.txt
#SBATCH -e /<path_to_working_directory>/step1/Log/logErrorAlignAll_ChIP.txt
#
# set the account for hpc cluster user
#SBATCH --account=userid
#
#SBATCH --export=ALL
#
########## BEGIN ACTUAL COMMANDS 


# Required modules
module purge
module load java/latest
module load samtools/1.2
module load python/2.7
module load cutadapt/1.8.1
module load picard-tools/1.92
module load bowtie/2.3.1

export BOWTIE2_INDEXES=/<path_to_bowtie2_index_directory>


# Begin commands for adapter trimming and alignment of Control H3K18ac ChIP files
# There are two repilicates per histone mark, and two pair-end PE fastq files per replicate
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H3K18ac_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H3K18ac_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Cnt1H3K18ac_PE1.fastq.gz /<path_to_working_directory>/fastq/Cnt1H3K18ac_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H3K18ac_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H3K18ac_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H3K18ac_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H3K18ac_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Cnt1H3K18ac_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Cnt1H3K18ac_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Cnt1H3K18ac_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Cnt1H3K18ac_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt1H3K18ac_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Cnt1H3K18ac_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Cnt1H3K18ac_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt1H3K18ac_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Cnt1H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H3K18ac_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H3K18ac_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Cnt1H3K18ac_PE1.fastq.gz /<path_to_working_directory>/fastq/Cnt1H3K18ac_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H3K18ac_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H3K18ac_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H3K18ac_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H3K18ac_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Cnt2H3K18ac_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Cnt2H3K18ac_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Cnt2H3K18ac_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Cnt2H3K18ac_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt2H3K18ac_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Cnt2H3K18ac_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Cnt2H3K18ac_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt2H3K18ac_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Cnt2H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam


#Begin commands for adapter trimming and alignment of Propionyl H3K18pr ChIP files
#There are two repilicates per histone mark, and two pair-end PE fastq files per replicate
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H3K18pr_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H3K18pr_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Prop1H3K18pr_PE1.fastq.gz /<path_to_working_directory>/fastq/Prop1H3K18pr_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H3K18pr_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H3K18pr_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H3K18pr_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H3K18pr_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Prop1H3K18pr_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Prop1H3K18pr_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Prop1H3K18pr_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Prop1H3K18pr_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Prop1H3K18pr_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Prop11H3K18ac_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Prop1H3K18pr_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Prop1H3K18pr_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Prop1H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H3K18pr_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H3K18pr_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Prop2H3K18pr_PE1.fastq.gz /<path_to_working_directory>/fastq/Prop1H3K18pr_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H3K18pr_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H3K18pr_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H3K18pr_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H3K18pr_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Prop1H3K18pr_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Prop1H3K18pr_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Prop2H3K18pr_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Prop2H3K18pr_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Prop2H3K18pr_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Prop2H3K18pr_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Prop1H3K18pr_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Prop2H3K18pr_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Prop2H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam


#Begin commands for adapter trimming and alignment of Control H4K12ac ChIP files
#There are two repilicates per histone mark, and two pair-end PE fastq files per replicate
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H4K12ac_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H4K12ac_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Cnt1H4K12ac_PE1.fastq.gz /<path_to_working_directory>/fastq/Cnt1H4K12ac_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H4K14ac_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H4K12ac_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H4K12ac_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt1H4K12ac_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Cnt1H4K12ac_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Cnt1H4K12ac_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Cnt1H3K18ac_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Cnt1H3K18ac_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt1H3K18ac_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Cnt1H4K12ac_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Cnt1H3K18ac_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt1H3K18ac_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Cnt1H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H3K18ac_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H3K18ac_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Cnt2H3K18ac_PE1.fastq.gz /<path_to_working_directory>/fastq/Cnt2H3K18ac_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H4K12ac_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H4K12ac_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H4K12ac_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Cnt2H4K12ac_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Cnt2H4K12ac_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Cnt2H3K18ac_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Cnt2H4K12ac_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Cnt2H4K12ac_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt2H4K12ac_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Cnt2H4K12ac_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Cnt2H3K18ac_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Cnt2H4K12ac_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Cnt2H4K12ac_sorted_dedup_filterUnmap_filterChrM_aln.bam


#Begin commands for adapter trimming and alignment of Propionyl H4K12pr ChIP files
#There are two repilicates per histone mark, and two pair-end PE fastq files per replicate
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H4K12pr_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H4K12pr_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Prop1H4K12pr_PE1.fastq.gz /<path_to_working_directory>/fastq/Prop1H4K12pr_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H4K12pr_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H4K12pr_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H4K12pr_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop1H4K12pr_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Prop1H4K12pr_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Prop1H4K12pr_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Prop1H4K12pr_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Prop1H4K12pr_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Prop1H4K12pr_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Prop1H4K12ac_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Prop1H4K12pr_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Prop1H4K12pr_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Prop1H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam

cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 30 -q 15 -o /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H4K12pr_trimmed_PE1.fastq.gz -p /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H4K12pr_trimmed_PE2.fastq.gz /<path_to_working_directory>/fastq/Prop2H4K12pr_PE1.fastq.gz /<path_to_working_directory>/fastq/Prop2H4K12pr_PE2.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H4K12pr_trimmed_PE1.fastq.gz
gunzip /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H4K12pr_trimmed_PE2.fastq.gz
bowtie2 -q --phred33 -X 2000 --fr -p 8 -x hg38 -1 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H4K12pr_trimmed_PE1.fastq  -2 /<path_to_working_directory>/step1/TrimmedFASTQ/Prop2H4K12pr_trimmed_PE2.fastq | samtools view -Sb - > /<path_to_working_directory>/step1/RawAlign/Prop2H4K12pr_raw_aln.bam
samtools sort -@ 8 -m 2G /<path_to_working_directory>/step1/RawAlign/Prop1H4K12pr_raw_aln.bam /<path_to_working_directory>/step1/SortedAlign/Prop2H4K12pr_sorted_aln
java -jar -Xms8g -Xmx8g /<path_to_picard-tools>/MarkDuplicates.jar INPUT=/<path_to_working_directory>/step1/SortedAlign/Prop2H4K12pr_sorted_aln.bam OUTPUT=/<path_to_working_directory>/step1/SortedAndDedupAlign/Prop2H4K12pr_sorted_dedup_aln.bam METRICS_FILE=/<path_to_working_directory>/step1/Temp/Prop2H4K12ac_temp_metrics.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true QUIET=true
rm /<path_to_working_directory>/step1/Temp/Prop2H4K12pr_temp_metrics.txt
samtools view -b -f 1 -F 12 -L /<path_to_working_directory>/inclusionZones_removeChrM.bed /<path_to_working_directory>/step1/SortedAndDedupAlign/Prop2H4K12pr_sorted_dedup_aln.bam > /<path_to_working_directory>/step1/FinalAlign/Prop2H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam


