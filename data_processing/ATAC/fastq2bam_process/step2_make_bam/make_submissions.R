input <- list.files("/group/wi38/data/ATACseq/181206/concat_files", pattern="R1", full.names=T)
submission_dir <- "/home/yuanh/analysis/WI38_ATAC/scripts/step2_make_bam/submissions/"
output_dir <- "/home/yuanh/analysis/WI38_ATAC/process/bam/"
dir.create(submission_dir, recursive = T)
dir.create(output_dir)

for (i in 1:length(input)) {
  myscript <- "/home/yuanh/programs/myscripts/NGS_scripts/fastq2bam_PE_bowtie.sh"
  index <- "/home/yuanh/programs/genomes/hg38/bowtie_index/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index"
  R1 <- input[i]
  R2 <- gsub("R1","R2",R1)
  adapter <- "/home/yuanh/programs/source/Trimmomatic-0.36/adapters/NexteraPE-PE.fa"
  sample <- paste0(output_dir, gsub(".*/|.fastq.gz","",R1))
  cores <- 1
  logfile <- paste0(output_dir, gsub(".*/|.fastq.gz","",R1), ".log")
  cmd <- paste(myscript, index, R1, R2, adapter, sample, cores, ">", logfile, "2>&1")
  
  submission_file <- paste0(submission_dir,  gsub(".*/|.fastq.gz","",R1), ".sh")
  
  # write slurm options
  write("#! /bin/bash", submission_file,append=FALSE)
  write("#SBATCH --ntasks=1", submission_file,append=TRUE)
  write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
  write("#SBATCH --mem=10G", submission_file,append=TRUE)
  write(paste0("#SBATCH --job-name=step2_make_bam_",i), submission_file,append=TRUE)
  write(paste0("#SBATCH --output=",gsub(".*/|.fastq.gz","",R1),".%j.out"), submission_file,append=TRUE)
  
  # write script
  write("module add bowtie2", submission_file,append=TRUE)
  write(cmd, submission_file,append=TRUE)
  
  # submit
  system(paste0("sbatch < ", submission_file))
}

