#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4.1_macs2pool_11
#SBATCH --output=PDL37_TP5_omni.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL37_TP5_omni_A_S36_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL37_TP5_omni_B_S37_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL37_TP5_omni_C_S38_SE.bed -f BED -n PDL37_TP5_omni --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
