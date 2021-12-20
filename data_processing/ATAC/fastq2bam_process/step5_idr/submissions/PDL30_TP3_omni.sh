#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_9
#SBATCH --output=PDL30_TP3_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL30_TP3_omni_A_S24_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL30_TP3_omni_B_S25_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL30_TP3_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL30_TP3_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL30_TP3_omni_A_S24_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL30_TP3_omni_C_S26_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL30_TP3_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL30_TP3_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL30_TP3_omni_B_S25_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL30_TP3_omni_C_S26_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL30_TP3_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL30_TP3_omni_BC_idr.txt
