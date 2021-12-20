#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name=step6_make_bw_34
#SBATCH --output=PDL25_TP2_std_B_S22.%j.out
Rscript /home/yuanh/programs/myscripts/NGS_scripts/makeBW_atac.R /home/yuanh/analysis/WI38_ATAC/process/bam/PDL25_TP2_std_B_S22_R1_001.filter.rmdup.bam 73
