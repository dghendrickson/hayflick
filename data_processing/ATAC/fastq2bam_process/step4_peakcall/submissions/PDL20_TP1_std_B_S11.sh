#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step4_macs2_28
#SBATCH --output=PDL20_TP1_std_B_S11.%j.out
source /home/yuanh/.bashrc
source activate py2
macs2 callpeak -t /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL20_TP1_std_B_S11_SE.bed -f BED -n PDL20_TP1_std_B_S11 --outdir /home/yuanh/analysis/WI38_ATAC/process/macs2out/ --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73
