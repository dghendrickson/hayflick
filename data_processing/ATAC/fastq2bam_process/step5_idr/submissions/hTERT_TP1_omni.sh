#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_1
#SBATCH --output=hTERT_TP1_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP1_omni_A_S1_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP1_omni_B_S2_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/hTERT_TP1_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/hTERT_TP1_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP1_omni_A_S1_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP1_omni_C_S3_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/hTERT_TP1_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/hTERT_TP1_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP1_omni_B_S2_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/hTERT_TP1_omni_C_S3_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/hTERT_TP1_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/hTERT_TP1_omni_BC_idr.txt
