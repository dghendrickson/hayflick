#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4_macs2_7
#SBATCH --output=hTERT_TP2_omni_A_S13.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP2_omni_A_S13_SE.bed -f BED -n hTERT_TP2_omni_A_S13 --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
