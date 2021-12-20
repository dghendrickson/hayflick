#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_38
#SBATCH --output=PDL30_TP3_omni_C_S26.%j.out
samtools sort -n -T PDL30_TP3_omni_C_S26_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL30_TP3_omni_C_S26_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL30_TP3_omni_C_S26_SE.bed
