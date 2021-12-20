#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8000
#SBATCH --job-name=step1
#SBATCH --output=step1_fastqc.%j.out

#notification options
#SBATCH --mail-type=end                        	# send email upon completion
#SBATCH --mail-user=yuanh@calicolabs.com       	# address to send status email

# script
# step1: run fastqc on each relevant fastq file
cd /home/yuanh/analysis/WI38_ATAC/process
mkdir raw
for ENTRY in /group/wi38/data/ATACseq/181206/concat_files/*_omni_*.gz
do
 fastqc -q -t 16 -o raw ${ENTRY}
 echo ${ENTRY}
done
