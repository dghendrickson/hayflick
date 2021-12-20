#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --job-name=step6_make_bw_32
#SBATCH --output=PDL25_TP2_omni_C_S18.%j.out
Rscript /home/yuanh/programs/myscripts/NGS_scripts/makeBW_atac.R /home/yuanh/analysis/WI38_ATAC/process/bam/PDL25_TP2_omni_C_S18_R1_001.filter.rmdup.bam 73
