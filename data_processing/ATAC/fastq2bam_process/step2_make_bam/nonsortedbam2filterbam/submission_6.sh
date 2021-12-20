#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=nonsort2bam_6
#SBATCH --output=PDL20_TP1_omni_C_S6_R1_001.%j.out
samtools sort -m 8G -o /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.filter.bam /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.nonSorted.bam
java -Xmx8g -jar /home/yuanh/programs/source/picard-2.20.3/picard.jar MarkDuplicates REMOVE_DUPLICATES=true AS=true I=/home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.filter.bam O=/home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.filter.rmdup.bam M=/home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001_rmdup_metric.txt
samtools index /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.filter.rmdup.bam
samtools flagstat /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.filter.rmdup.bam > /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.flagtat
rm /home/yuanh/analysis/WI38_ATAC/process/bam/PDL20_TP1_omni_C_S6_R1_001.filter.bam
