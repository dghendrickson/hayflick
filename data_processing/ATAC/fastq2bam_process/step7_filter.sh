#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8000
#SBATCH --job-name=step7
#SBATCH --output=step7_filter.%j.out

#notification options
#SBATCH --mail-type=end                         # send email upon completion
#SBATCH --mail-user=yuanh@calicolabs.com        # address to send status email

# script
# step7: filter idr files
for ENTRY in /home/yuanh/analysis/WI38_ATAC/process/idr/*.txt
do
 Rscript /home/yuanh/programs/myscripts/NGS_scripts/filter.R ${ENTRY} hg38 0.05
 echo ${ENTRY}
done
