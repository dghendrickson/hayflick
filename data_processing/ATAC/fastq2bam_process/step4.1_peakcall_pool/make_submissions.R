######################
# macs2 on pool data #
######################
input <- list.files("/home/yuanh/analysis/WI38_ATAC/process/macs2out", pattern=".narrowPeak", full.names=T)
input <- input[grep("_omni_", input)]
input <- unique(gsub(".*/|_[A|B|C]_.*","",input))

submission_dir <- "/home/yuanh/analysis/WI38_ATAC/scripts/step4.1_peakcall_pool/submissions/"
output_dir <- "/home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/"
dir.create(submission_dir, recursive = T)
dir.create(output_dir)

for (i in 1:length(input)) {
  
  replicates <- list.files("/home/yuanh/analysis/WI38_ATAC/process/shift_bed", pattern=input[i], full.names=T)
  
  cmd <- paste0("macs2 callpeak -t ", 
                paste(replicates, collapse = " "), 
                " -f BED -n ",input[i],
                " --outdir ",output_dir,
                " --keep-dup all -g hs -p 1e-1 --nomodel --shift -37 --extsize 73")
  
  submission_file <- paste0(submission_dir,  input[i], ".sh")
  
  # write slurm options
  write("#! /bin/bash", submission_file,append=FALSE)
  write("#SBATCH --ntasks=1", submission_file,append=TRUE)
  write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
  write("#SBATCH --mem=10G", submission_file,append=TRUE)
  write(paste0("#SBATCH --job-name=step4.1_macs2pool_",i), submission_file,append=TRUE)
  write(paste0("#SBATCH --output=",input[i],".%j.out"), submission_file,append=TRUE)
  
  # write script
  write("source /home/yuanh/.bashrc", submission_file, append=TRUE)
  write("source activate py2", submission_file, append=TRUE)
  write(cmd, submission_file,append=TRUE)
  
  # submit
  system(paste0("sbatch < ", submission_file))
}

