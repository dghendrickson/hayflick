#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_7
#SBATCH --output=PDL20_TP1_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL20_TP1_omni_A_S4_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL20_TP1_omni_B_S5_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL20_TP1_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL20_TP1_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL20_TP1_omni_A_S4_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL20_TP1_omni_C_S6_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL20_TP1_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL20_TP1_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL20_TP1_omni_B_S5_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL20_TP1_omni_C_S6_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL20_TP1_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL20_TP1_omni_BC_idr.txt
