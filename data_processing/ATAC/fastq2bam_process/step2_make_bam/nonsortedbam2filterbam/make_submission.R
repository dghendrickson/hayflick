# 11 large sequencing files failed after alignment, at the rm .fq step.
# the exact reason for the failure is not clear. The error is "stale file handle".
# Here is the script to continue with sort, rmdup.

input <- list.files("/home/yuanh/analysis/WI38_ATAC/process/bam", pattern="nonSorted.bam", full.names=T)
submission_dir <- "/home/yuanh/analysis/WI38_ATAC/scripts/step2_make_bam/nonsortedbam2filterbam/"

for (i in 1:length(input)) {
  submission_file <- paste0(submission_dir, "submission_", i, ".sh")
  
  # write slurm options
  write("#! /bin/bash", submission_file,append=FALSE)
  write("#SBATCH --ntasks=1", submission_file,append=TRUE)
  write("#SBATCH --cpus-per-task=1", submission_file,append=TRUE)
  write("#SBATCH --mem=10G", submission_file,append=TRUE)
  write(paste0("#SBATCH --job-name=nonsort2bam_",i), submission_file,append=TRUE)
  write(paste0("#SBATCH --output=",gsub(".*/|.nonSorted.bam","",input[i]),".%j.out"), submission_file,append=TRUE)
  
  write(paste("samtools sort -m 8G -o", gsub(".nonSorted.bam",".filter.bam",input[i]), input[i]),
        submission_file,append=TRUE)
  write(paste0("java -Xmx8g -jar /home/yuanh/programs/source/picard-2.20.3/picard.jar MarkDuplicates REMOVE_DUPLICATES=true AS=true I=",
         gsub(".nonSorted.bam",".filter.bam",input[i]),
         " O=", gsub(".nonSorted.bam",".filter.rmdup.bam",input[i]),
         " M=", gsub(".nonSorted.bam","_rmdup_metric.txt",input[i])),
        submission_file,append=TRUE)
        
  write(paste("samtools index", gsub(".nonSorted.bam",".filter.rmdup.bam",input[i])),
        submission_file,append=TRUE)
  write(paste("samtools flagstat",gsub(".nonSorted.bam",".filter.rmdup.bam",input[i]),">",gsub(".nonSorted.bam",".flagtat",input[i])),
        submission_file,append=TRUE)
  write(paste("rm",  gsub(".nonSorted.bam",".filter.bam",input[i])),
        submission_file,append=TRUE)
  
  # submit
  system(paste0("sbatch < ", submission_file))
}

