#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_44
#SBATCH --output=PDL37_TP5_omni_C_S38.%j.out
samtools sort -n -T PDL37_TP5_omni_C_S38_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL37_TP5_omni_C_S38_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL37_TP5_omni_C_S38_SE.bed
