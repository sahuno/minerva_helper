#!/bin/bash

#date: March 30th  2021
#make read groups from demultiplexed fastq files 

# To run use this program, 
# $ source get_read_group_fastq.sh <directory/to/fastq/files/> -O cohort_or_study_name


echo "Using $1 as fastq directory"
echo "Using $2 as project name";

##create tmp variables;
read_groups_all=$2

#create file to be updated with each samples read group in, PU names
echo -e "Sample_name" ' \t ' "Read_R1_or_R2" ' \t ' "Platform_Unit_PU" ' \t ' "Read_Group_ID_RGID" ' \t ' "Library_LB" ' \t ' "Lane" ' \t ' "File_path_on_drive" ' \t ' "Index_seq" ' \t ' "file_Name">> ${read_groups_all}.txt


##this is how a fastq file header looks like
general_Fastq_file_header="$(printf '%s' '@instrument_id:run_number:flowcell_ID:lane:tile:x_pos:y_pos member_read_pair:is_filtered:control_number:index_sequence')"
#FLOWCELL_BARCODE = @(instrument id):(run number):(flowcell ID)
IFS=': ' read -a general_Fastq_file_header_array <<< $general_Fastq_file_header
printf "%s\t" ${general_Fastq_file_header_array[@]}

## see these references
#https://angus.readthedocs.io/en/2017/Read_group_info.html
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups



for file in ${1}*.fastq.gz
do
echo "Processing $file"
#tmp_header="tmp_header"
tmp_header=$(zcat $file | head -n 1) #save header into temp vriable ($tmp_header); should look something like this $temp_header=$(zcat $file | head -n 1)

#convert fastq header into array
IFS=': ' read -a header_array <<< $tmp_header
echo ${header_array[@]} # do this to see contents of header; echo ${header_array[0]} is instrument ID, echo ${header_array[10]} is the index sequence

## delete "@" character of 
delete=$(echo "@")
header_array=("${header_array[@]/$delete}")

#begin assigning variables, follow this format
file_Path=$file
file_Name_toAppend=$(basename $file)
SM=$(basename -s .fastq.gz $file | cut -d '_' -f 1) ##please remember to trim away any "R1/R2" to match sample name on terra tables, I assume file names are sperated by delimiter '_'
lane=${header_array[3]}
LB=$SM.${header_array[10]}
read=${header_array[7]}
flowcell_barcode=${header_array[0]}.${header_array[1]}.${header_array[2]}
PU=$flowcell_barcode.${header_array[3]}.${header_array[10]}
RGID=$flowcell_barcode.${header_array[3]}
#append extracted information to file
echo -e "$SM" ' \t ' "$read" ' \t ' "$PU" ' \t ' "$RGID" ' \t ' "$LB" ' \t ' "$lane" ' \t ' "$file_Path" ' \t ' "${header_array[10]}" ' \t ' "$file_Name_toAppend">> ${read_groups_all}.txt
done


