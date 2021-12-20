#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_11
#SBATCH --output=PDL37_TP5_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL37_TP5_omni_A_S36_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL37_TP5_omni_B_S37_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL37_TP5_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL37_TP5_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL37_TP5_omni_A_S36_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL37_TP5_omni_C_S38_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL37_TP5_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL37_TP5_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL37_TP5_omni_B_S37_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL37_TP5_omni_C_S38_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL37_TP5_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL37_TP5_omni_BC_idr.txt
