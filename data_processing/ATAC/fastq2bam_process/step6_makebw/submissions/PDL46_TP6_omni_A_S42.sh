#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name=step6_make_bw_45
#SBATCH --output=PDL46_TP6_omni_A_S42.%j.out
Rscript /home/yuanh/programs/myscripts/NGS_scripts/makeBW_atac.R /home/yuanh/analysis/WI38_ATAC/process/bam/PDL46_TP6_omni_A_S42_R1_001.filter.rmdup.bam 73
