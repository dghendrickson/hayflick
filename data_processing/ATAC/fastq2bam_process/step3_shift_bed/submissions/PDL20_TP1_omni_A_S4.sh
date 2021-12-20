#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_24
#SBATCH --output=PDL20_TP1_omni_A_S4.%j.out
samtools sort -n -T PDL20_TP1_omni_A_S4_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_A_S4_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL20_TP1_omni_A_S4_SE.bed
