#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_42
#SBATCH --output=PDL37_TP5_omni_A_S36.%j.out
samtools sort -n -T PDL37_TP5_omni_A_S36_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL37_TP5_omni_A_S36_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL37_TP5_omni_A_S36_SE.bed
