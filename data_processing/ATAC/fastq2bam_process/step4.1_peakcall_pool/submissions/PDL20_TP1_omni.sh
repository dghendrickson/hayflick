#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4.1_macs2pool_7
#SBATCH --output=PDL20_TP1_omni.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL20_TP1_omni_A_S4_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL20_TP1_omni_B_S5_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL20_TP1_omni_C_S6_SE.bed -f BED -n PDL20_TP1_omni --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
