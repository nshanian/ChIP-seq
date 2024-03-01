#!/bin/bash -l
#
# Set the name of the job
#SBATCH --job-name=chip_dbind
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
#SBATCH -n 1
#
# Set output and error log files
#SBATCH -o /<path_to_working_directory>/diffbind/out.txt
#SBATCH -e /<path_to_working_directory>/diffbind/error.txt
#
# set the account for hpc cluster user
#SBATCH --account=userid
#
#
#SBATCH --export=ALL


# set diffbind directory
dir=/<path_to_working_directory>/diffbind
sampleSheet=<path_to_working_directory>/diffbind/H4K12pr_diffbind.csv
name=chip

########## BEGIN ACTUAL COMMANDS

#Required modules
module load r/3.6



#Begin commands
Rscript <path_to_working_directory>/diffbind/diffbind_chip.R ${dir} ${sampleSheet} ${name}

