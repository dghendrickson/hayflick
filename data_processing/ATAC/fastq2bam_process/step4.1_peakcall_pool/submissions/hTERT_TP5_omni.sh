#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4.1_macs2pool_4
#SBATCH --output=hTERT_TP5_omni.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP5_omni_A_S33_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP5_omni_B_S34_SE.bed /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP5_omni_C_S35_SE.bed -f BED -n hTERT_TP5_omni --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
