input <- list.files("/home/yuanh/analysis/WI38_ATAC/process/bam", pattern=".filter.rmdup.bam$", full.names=T)
submission_dir <- "/home/yuanh/analysis/WI38_ATAC/scripts/step6_makebw/submissions/"
dir.create(submission_dir, recursive = T)

for (i in 1:length(input)) {
  cmd  <- paste0("Rscript /home/yuanh/programs/myscripts/NGS_scripts/makeBW_atac.R ", input[i], " 73")
  
  # write slurm options
  submission_file <- paste0(submission_dir,  gsub(".*/|_R1_001.filter.rmdup.bam","",input[i]), ".sh")
  write("#! /bin/bash", submission_file,append=FALSE)
  write("#SBATCH --ntasks=1", submission_file,append=TRUE)
  write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
  write("#SBATCH --mem=40G", submission_file,append=TRUE)
  write(paste0("#SBATCH --job-name=step6_make_bw_",i), submission_file,append=TRUE)
  write(paste0("#SBATCH --output=", gsub(".*/|_R1_001.filter.rmdup.bam","",input[i]),".%j.out"), submission_file,append=TRUE)
  
  # write script
  write(cmd, submission_file,append=TRUE)
  
  # submit
  system(paste0("sbatch < ", submission_file))
}

