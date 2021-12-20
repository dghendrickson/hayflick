#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4.1_macs2pool_13
#SBATCH --output=PDL50_TP7_omni.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL50_TP7_omni_A_S48_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL50_TP7_omni_B_S49_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL50_TP7_omni_C_S50_SE.bed -f BED -n PDL50_TP7_omni --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
