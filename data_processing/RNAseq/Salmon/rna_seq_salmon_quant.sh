#!/bin/sh

##Salmon vesion 0.82
## used cellRanger 3.0 hg38 reference and fasta see cellRanger_make_reference_and_fasta_.sh
## Used Salmon "index" to generate Salmon index with default settings

tx_index=/path/to/transcript/reference


for fastq in ls .
do
fast2=${fastq/%R1*/R2_001.fastq.gz}
       srun -p standard   --mem=60000 -e ./salmon.err  salmon quant -i $tx_index -l ISR -1 $fastq   -2  $fast2 -o  ./${fastq%.*}.squant --useVBOpt --numBootstraps 30 &

done

