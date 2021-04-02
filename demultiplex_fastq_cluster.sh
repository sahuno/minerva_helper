#!/bin/bash

#!/bin/bash
#BSUB -J  %J  # Job name
#BSUB -P acc_polakp02a # allocation account
#BSUB -q long # queue
#BSUB -n 4 # number of compute cores
#BSUB -W 0:30 # walltime in HH:MM
#BSUB -R rusage[mem=12000] # 12 GB of memory (4 GB per core)
#BSUB -oo stdout # output log (%J : JobID)
#BSUB -eo stderr # error log
#BSUB -L /bin/bash # Initialize the execution environment

## how to run is give it path containg fastq and they will be demultiplex, compress after wards and aling with approaproatet read groups
#fastq_dir=$1 #use this interactive mode
fastq_dir=/sc/arion/projects/canan/data_ahuno/CRC_PR_03192021/fastq/DrMarciaCruz_WES_CRC_LateOnset/Tumor/
#line=/sc/arion/projects/canan/data_ahuno/CRC_PR_03192021/fastq/DrMarciaCruz_WES_CRC_EarlyOnset/Tumor/test_pipelines/SC173887_AACTCACC_L001_R1_001_HQ_paired.fastq.gz



for line in ${fastq_dir}*.fastq.gz
do

echo "processing fastq $line "
Sample_name=$(basename -s .fastq.gz "$line")

time zcat "$line" |
        awk -v sname="$Sample_name" '
            BEGIN {FS = ":"} 
            {
                lane=$4;
                fileName = sname".lane."lane".fastq" 
                print > fileName
                for (i = 1; i <= 3; i++) {
                    getline
                    print > fileName
                    
                }
            }'
done


#for 
# real    1m43.400s
# user    2m40.082s
# sys     0m17.742s

#for y in ${fastq_dir}*.fastq
#Sample_name2Comp=$(basename -s .fastq "$y")
#do
#echo "compressing $y"
#time gzip $y > ${Sample_name2Comp}.gz
#done
