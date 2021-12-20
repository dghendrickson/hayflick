#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_40
#SBATCH --output=PDL33_TP4_omni_B_S31.%j.out
samtools sort -n -T PDL33_TP4_omni_B_S31_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL33_TP4_omni_B_S31_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL33_TP4_omni_B_S31_SE.bed
