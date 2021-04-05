#!/bin/bash
#BSUB -J  compression  # Job name
#BSUB -P acc_polakp02a # allocation account
#BSUB -q long # queue
#BSUB -n 6 # number of compute cores
#BSUB -W 1:00 # walltime in HH:MM
#BSUB -R rusage[mem=12000] # 12 GB of memory (4 GB per core)
#BSUB -oo stdout # output log (%J : JobID)
#BSUB -eo stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

##date, March 31st, to compress fastqs


while read line
  do
      echo "$line"
      Sample_name2Comp=$(basename -s .fastq "$line")
      echo "compressing $line"
      gzip $line
done < /sc/arion/projects/canan/data_ahuno/CRC_PR_03192021/fastq/DrMarciaCruz_WES_CRC_EarlyOnset/Tumor/demuliplexed_fastq/list_fastq_files_1_comp.txt
