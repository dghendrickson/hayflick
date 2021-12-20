#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_45
#SBATCH --output=PDL46_TP6_omni_A_S42.%j.out
samtools sort -n -T PDL46_TP6_omni_A_S42_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL46_TP6_omni_A_S42_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL46_TP6_omni_A_S42_SE.bed
