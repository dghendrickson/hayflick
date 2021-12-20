#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_27
#SBATCH --output=PDL20_TP1_std_A_S10.%j.out
samtools sort -n -T PDL20_TP1_std_A_S10_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_std_A_S10_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL20_TP1_std_A_S10_SE.bed
