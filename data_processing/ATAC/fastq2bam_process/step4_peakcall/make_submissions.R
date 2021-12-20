input <- list.files("/home/yuanh/analysis/WI38_ATAC/process/shift_bed", pattern="_SE.bed", full.names=T)
submission_dir <- "/home/yuanh/analysis/WI38_ATAC/scripts/step4_peakcall/submissions/"
output_dir <- "/home/yuanh/analysis/WI38_ATAC/process/macs2out/"
dir.create(submission_dir, recursive = T)
dir.create(output_dir)

for (i in 1:length(input)) {
  sample_id <- gsub(".*/|_SE.bed","",input[i])
  
  
  cmd <- paste0("macs2 callpeak -t ",input[i], 
                " -f BED -n ",sample_id,
                " --outdir ",output_dir,
                " --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73")
  
  
  submission_file <- paste0(submission_dir,  sample_id, ".sh")
  
  # write slurm options
  write("#! /bin/bash", submission_file,append=FALSE)
  write("#SBATCH --ntasks=1", submission_file,append=TRUE)
  write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
  write("#SBATCH --mem=10G", submission_file,append=TRUE)
  write(paste0("#SBATCH --job-name=step4_macs2_",i), submission_file,append=TRUE)
  write(paste0("#SBATCH --output=",sample_id,".%j.out"), submission_file,append=TRUE)
  
  # write script
  write("source /home/yuanh/.bashrc", submission_file, append=TRUE)
  write("source activate py2", submission_file, append=TRUE)
  write(cmd, submission_file,append=TRUE)
  
  # submit
  system(paste0("sbatch < ", submission_file))
}

