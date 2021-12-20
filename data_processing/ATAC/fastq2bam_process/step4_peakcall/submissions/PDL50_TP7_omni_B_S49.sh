#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4_macs2_49
#SBATCH --output=PDL50_TP7_omni_B_S49.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL50_TP7_omni_B_S49_SE.bed -f BED -n PDL50_TP7_omni_B_S49 --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
