#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_2
#SBATCH --output=hTERT_TP2_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP2_omni_A_S13_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP2_omni_B_S14_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/hTERT_TP2_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/hTERT_TP2_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP2_omni_A_S13_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP2_omni_C_S15_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/hTERT_TP2_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/hTERT_TP2_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP2_omni_B_S14_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP2_omni_C_S15_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/hTERT_TP2_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/hTERT_TP2_omni_BC_idr.txt
