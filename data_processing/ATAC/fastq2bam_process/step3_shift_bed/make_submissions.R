input <- list.files("/home/yuanh/analysis/WI38_ATAC/process/bam", pattern=".filter.rmdup.bam$", full.names=T)
submission_dir <- "/home/yuanh/analysis/WI38_ATAC/scripts/step3_shift_bed/submissions/"
output_dir <- "/home/yuanh/analysis/WI38_ATAC/process/shift_bed/"
dir.create(submission_dir, recursive = T)
dir.create(output_dir)

for (i in 1:length(input)) {
  sample_id <- gsub(".*/|.filter.rmdup.bam|_R1_001","",input[i])
  
  cmd <- paste0("samtools sort -n -T ",sample_id,"_nsorted ", input[i],
                " |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > ",
                output_dir,sample_id,"_SE.bed")
  
  submission_file <- paste0(submission_dir,  sample_id, ".sh")
  
  # write slurm options
  write("#! /bin/bash", submission_file,append=FALSE)
  write("#SBATCH --ntasks=1", submission_file,append=TRUE)
  write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
  write("#SBATCH --mem=10G", submission_file,append=TRUE)
  write(paste0("#SBATCH --job-name=step3_shift_bed_",i), submission_file,append=TRUE)
  write(paste0("#SBATCH --output=",sample_id,".%j.out"), submission_file,append=TRUE)
  
  # write script
  write(cmd, submission_file,append=TRUE)
  
  # submit
  system(paste0("sbatch < ", submission_file))
}

