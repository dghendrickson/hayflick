#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_13
#SBATCH --output=hTERT_TP4_omni_B_S28.%j.out
samtools sort -n -T hTERT_TP4_omni_B_S28_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/hTERT_TP4_omni_B_S28_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP4_omni_B_S28_SE.bed
