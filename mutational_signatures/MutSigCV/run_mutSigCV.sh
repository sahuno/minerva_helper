#!/bin/bash

##MutSigCV; March 10th 2021 by Samuel Ahuno
#Author; Samuel Ahuno

#BSUB -J "mutsig_projectA" ##job name
#BSUB -P project_Acc_Name
#BSUB -q premium
#BSUB -n 2
#BSUB -W 1:00 # I like running suing short times, cos that gets me ahead of the queue
#BSUB -R rusage[mem=4000]
#BSUB -R himem
#BSUB -o %J.stdout
#BSUB -eo %J.stderr
#BSUB -L /bin/bash

module load mutsig matlab/R2016b

MutSigCV /path/to/maf_file_in_hg19_format.maf \
/hpc/users/ahunos01/apps/MutSigCV_1.41/exdata/exome_full192.coverage.txt \
/hpc/users/ahunos01/apps/MutSigCV_1.41/exdata/gene.covariates.txt \
basename_of_output_file /hpc/users/ahunos01/apps/MutSigCV_1.41/exdata/mutation_type_dictionary_file.txt \
/hpc/users/ahunos01/apps/MutSigCV_1.41/exdata/chr_files_hg19

