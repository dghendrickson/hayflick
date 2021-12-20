#######
# IDR #
#######
input <- list.files("/home/yuanh/analysis/WI38_ATAC/process/macs2out", pattern=".narrowPeak", full.names=T)
input <- input[grep("_omni_", input)]
input <- unique(gsub(".*/|_[A|B|C]_.*","",input))

submission_dir <- "/home/yuanh/analysis/WI38_ATAC/scripts/step5_idr/submissions/"
output_dir <- "/home/yuanh/analysis/WI38_ATAC/process/idr/"
dir.create(submission_dir, recursive = T)
dir.create(output_dir)

for (i in 1:length(input)) {
  replicates <- list.files("/home/yuanh/analysis/WI38_ATAC/process/macs2out", pattern=input[i], full.names=T)
  replicates <- replicates[grep("narrowPeak", replicates)]
  
  pairwise <- combn(replicates, 2)
  
  cmds <- apply(pairwise, 2, function(x) {
    samples <- paste(x, collapse=" ")
    ids <- paste(gsub("_.*","",gsub(".*_omni_","",x)), collapse="")
    cmd <- paste0("idr --samples ", samples, 
                  " --peak-list ", "/home/yuanh/analysis/WI38_ATAC/process/macs2out_pool/", input[i], "_peaks.narrowPeak",
                  " --plot ",
                  " --output-file ",output_dir, input[i], "_", ids, "_idr.txt")
  })
  
  submission_file <- paste0(submission_dir,  input[i], ".sh")
  
  # write slurm options
  write("#! /bin/bash", submission_file,append=FALSE)
  write("#SBATCH --ntasks=1", submission_file,append=TRUE)
  write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
  write("#SBATCH --mem=10G", submission_file,append=TRUE)
  write(paste0("#SBATCH --job-name=step5_idr_",i), submission_file,append=TRUE)
  write(paste0("#SBATCH --output=",input[i],".%j.out"), submission_file,append=TRUE)
  
  # write script
  for (j in 1:length(cmds)) {
    write(cmds[j], submission_file,append=TRUE)
  }
  
  # submit
  system(paste0("sbatch < ", submission_file))
}

