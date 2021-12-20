#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4.1_macs2pool_6
#SBATCH --output=hTERT_TP7_omni.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP7_omni_A_S45_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP7_omni_B_S46_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP7_omni_C_S47_SE.bed -f BED -n hTERT_TP7_omni --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
