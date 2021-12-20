#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step5_idr_12
#SBATCH --output=PDL46_TP6_omni.%j.out
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL46_TP6_omni_A_S42_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL46_TP6_omni_B_S43_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL46_TP6_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL46_TP6_omni_AB_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL46_TP6_omni_A_S42_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL46_TP6_omni_C_S44_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL46_TP6_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL46_TP6_omni_AC_idr.txt
idr --samples /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL46_TP6_omni_B_S43_peaks.narrowPeak /home/yuanh/analysis/WI38_ATAC/process/macs2out/PDL46_TP6_omni_C_S44_peaks.narrowPeak --peak-list /home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/PDL46_TP6_omni_peaks.narrowPeak --plot  --output-file /home/yuanh/analysis/WI38_ATAC/process/idr/PDL46_TP6_omni_BC_idr.txt
