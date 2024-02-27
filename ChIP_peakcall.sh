#!/bin/bash -l
#
# Set the name of the job
#SBATCH --job-name=ChIPpeakcall
#
# Set the maximum memory allowed
#SBATCH --mem=30G
#
# Set the maximum run time
#SBATCH -t 48:00:00
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email
#
# The number of threads we will require
#SBATCH -n 1
#
# Set output and error log files
#SBATCH -o /<path_to_working_directory>/step2/Log/logOutPeakCall_ChIP.txt
#SBATCH -e /<path_to_working_directory>/step2/Log/logErrorPeakCall_ChIP.txt
#
# set the account for hpc cluster user
#SBATCH --account=userid
#
#SBATCH --export=ALL
#
########## BEGIN ACTUAL COMMANDS 


# Required modules
module load python/2.7
module load MACS2/2.1.0
module load samtools/1.2
module load java/latest
module load ucsc_tools/3.0.9


# Set necessary environmental variables for align2rawsignal
MCRROOT=/home/user/MCR/v714
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/runtime/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64
MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}
XAPPLRESDIR=${MCRROOT}/X11/app-defaults
export LD_LIBRARY_PATH
export XAPPLRESDIR
export MCR_CACHE_ROOT=/<path_to_working_directory>/align2rawsignal/temp


# Begin commands for peak calling and signal track generation for Control H3K18ac ChIP aligned files
macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Cnt1H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Cnt1H3K18ac -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt1H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H3K18ac.bedGraph" -of="bg" -n=5 -l=195 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H3K18ac.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H3K18ac.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H3K18ac.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt1H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Cnt1H3K18ac.mat" -of="mat" -n=5 -l=195 -k=epanechnikov -w=150 -f=0

macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Cnt2H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Cnt2H3K18ac -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt2H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H3K18ac.bedGraph" -of="bg" -n=5 -l=175 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H3K18ac.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H3K18ac.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H3K18ac.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt2H3K18ac_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Cnt2H3K18ac.mat" -of="mat" -n=5 -l=175 -k=epanechnikov -w=150 -f=0

# Begin commands for peak calling and signal track generation for Propionyl H3K18pr ChIP aligned files
macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Prop1H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Prop1H3K18pr -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop1H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H3K18pr.bedGraph" -of="bg" -n=5 -l=170 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H3K18pr.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H3K18pr.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H3K18pr.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop1H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Prop1H3K18pr.mat" -of="mat" -n=5 -l=170 -k=epanechnikov -w=150 -f=0

macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Prop2H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Prop2H3K18pr -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop2H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H3K18pr.bedGraph" -of="bg" -n=5 -l=160 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H3K18pr.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H3K18pr.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H3K18pr.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop2H3K18pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Prop2H3K18pr.mat" -of="mat" -n=5 -l=160 -k=epanechnikov -w=150 -f=0

# Begin commands for peak calling and signal track generation for Control H4K12ac ChIP aligned files
macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Cnt1H4K12ac_sorted_dedup_filterUnmap_filterChrM_aln -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Cnt1H4K12ac -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt1H4K12ac_sorted_dedup_filterUnmap_filterChrM_aln" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H4K12ac.bedGraph" -of="bg" -n=5 -l=160 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H4K12ac.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H4K12ac.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt1H4K12ac.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt1H4K12ac_sorted_dedup_filterUnmap_filterChrM_aln" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Cnt1H4K12ac.mat" -of="mat" -n=5 -l=160 -k=epanechnikov -w=150 -f=0

macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Cnt2H4K12ac_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Cnt2H4K12ac -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt2H4K12ac_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H4K12ac.bedGraph" -of="bg" -n=5 -l=140 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H4K12ac.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H4K12ac.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Cnt2H4K12ac.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Cnt2H4K12ac_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Cnt2H4K12ac.mat" -of="mat" -n=5 -l=140 -k=epanechnikov -w=150 -f=0

# Begin commands for peak calling and signal track generation for Propionyl H4K12pr ChIP aligned files
macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Prop1H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Prop1H4K12pr -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop1H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H4K12pr.bedGraph" -of="bg" -n=5 -l=170 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H4K12pr.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H4K12pr.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Prop1H4K12pr.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop1H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Prop1H4K12pr.mat" -of="mat" -n=5 -l=170 -k=epanechnikov -w=150 -f=0

macs2 callpeak -t /<path_to_working_directory>/step1/FinalAlign/Prop2H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam -g hs -f BAM --outdir /<path_to_working_directory>/step2/PeakFiles/ -n Prop2H4K12pr -q 0.05 --nomodel --extsize 75
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop2H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H4K12pr.bedGraph" -of="bg" -n=5 -l=225 -k=epanechnikov -w=150 -f=0
bedGraphToBigWig /<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H4K12pr.bedGraph /<path_to_working_directory>/hg38/hg38.genomeSizes_subset /<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H4K12pr.bw
rm /<path_to_working_directory>/step2/SignalTrackBigWig/Prop2H4K12pr.bedGraph
/<path_to_working_directory>/scripts/align2rawsignal/bin/align2rawsignal -i="/<path_to_working_directory>/step1/FinalAlign/Prop2H4K12pr_sorted_dedup_filterUnmap_filterChrM_aln.bam" -s="/<path_to_working_directory>/hg38/chr_fasta" -u="/<path_to_working_directory>/hg38/mappability/globalmap_k20tok101" -o="/<path_to_working_directory>/step2/SignalTrackMAT/Prop2H4K12pr.mat" -of="mat" -n=5 -l=225 -k=epanechnikov -w=150 -f=0


