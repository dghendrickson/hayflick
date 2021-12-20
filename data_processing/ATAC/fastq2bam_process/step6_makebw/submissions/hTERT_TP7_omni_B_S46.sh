#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name=step6_make_bw_22
#SBATCH --output=hTERT_TP7_omni_B_S46.%j.out
Rscript /home/yuanh/programs/myscripts/NGS_scripts/makeBW_atac.R /home/yuanh/analysis/WI38_ATAC/process/bam/hTERT_TP7_omni_B_S46_R1_001.filter.rmdup.bam 73
