#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step2_make_bam_19
#SBATCH --output=hTERT_TP6_omni_B_S40_R1_001.%j.out
module add bowtie2
/home/yuanh/programs/myscripts/NGS_scripts/fastq2bam_PE_bowtie.sh /home/yuanh/programs/genomes/hg38/bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index /group/wi38/data/ATACseq/181206/concat_files/hTERT_TP6_omni_B_S40_R1_001.fastq.gz /group/wi38/data/ATACseq/181206/concat_files/hTERT_TP6_omni_B_S40_R2_001.fastq.gz /home/yuanh/programs/source/Trimmomatic-0.36/adapters/NexteraPE-PE.fa /home/yuanh/analysis/WI38_ATAC/process/bam/hTERT_TP6_omni_B_S40_R1_001 1 > /home/yuanh/analysis/WI38_ATAC/process/bam/hTERT_TP6_omni_B_S40_R1_001.log 2>&1
