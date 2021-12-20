#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_8
#SBATCH --output=PDL25_TP2_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL25_TP2_omni_A_S16_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL25_TP2_omni_B_S17_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL25_TP2_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL25_TP2_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL25_TP2_omni_A_S16_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL25_TP2_omni_C_S18_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL25_TP2_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL25_TP2_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL25_TP2_omni_B_S17_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL25_TP2_omni_C_S18_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL25_TP2_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL25_TP2_omni_BC_idr.txt
