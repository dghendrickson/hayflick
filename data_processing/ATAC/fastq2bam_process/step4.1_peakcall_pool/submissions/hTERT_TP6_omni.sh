#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4.1_macs2pool_5
#SBATCH --output=hTERT_TP6_omni.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP6_omni_A_S39_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP6_omni_B_S40_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP6_omni_C_S41_SE.bed -f BED -n hTERT_TP6_omni --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
