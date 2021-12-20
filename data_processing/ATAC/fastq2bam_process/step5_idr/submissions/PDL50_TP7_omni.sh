#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_13
#SBATCH --output=PDL50_TP7_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL50_TP7_omni_A_S48_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL50_TP7_omni_B_S49_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL50_TP7_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL50_TP7_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL50_TP7_omni_A_S48_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL50_TP7_omni_C_S50_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL50_TP7_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL50_TP7_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL50_TP7_omni_B_S49_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL50_TP7_omni_C_S50_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL50_TP7_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL50_TP7_omni_BC_idr.txt
