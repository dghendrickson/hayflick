#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4.1_macs2pool_1
#SBATCH --output=hTERT_TP1_omni.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP1_omni_A_S1_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP1_omni_B_S2_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP1_omni_C_S3_SE.bed -f BED -n hTERT_TP1_omni --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
