#!/bin/bash

#date: March 30th  2021
#make read groups from fastq files 

## to run use this program, 
# $ source test.create.read.grp.cluster.sh -I /sc/arion/projects/canan/data_ahuno/CRC_PR_03192021/fastq/DrMarciaCruz_WES_CRC_LateOnset/Tumor/ -O cRC_test_read_groups

#while getopts I:P: flag
#do
#    case "${flag}" in
#        I) file_path=${OPTARG};;
#        P) project_name=${OPTARG};;
#    esac
#done

#echo "The file path is: $file_path";
#echo "The project name: $project_name";

##create tmp variables;
read_groups_all=$2
#create file to be updated with each samples read group in, PU names
echo -e " Sample_name" ' \t ' "Read_R1_or_R2" ' \t ' "Platform_Unit_PU" ' \t ' "Read_Group_ID_RGID" ' \t ' "Library_LB" ' \t ' "Lane" ' \t ' "File_path_on_drive" ' \t ' "Index_seq">> ${read_groups_all}.txt


#FILES=/Volumes/CRC_PR/DrMarciaCruz_WES_CRC_LateOnset/Tumor/
for file in ${1}*.fastq.gz
do
#tmp_header="tmp_header"
echo "Processing $file"


tmp_header=$(zcat $file | head -n 1) #save header into temp vriable ($tmp_header); should look something like this $temp_header=$(zcat $file | head -n 1)

#convert fastq header into array
IFS=': ' read -a header_array <<< $tmp_header
echo ${header_array[@]} # do this to see contents of header; echo ${header_array[0]} is instrument ID, echo ${header_array[10]} is the index sequence


## delete "@" character of 
delete=$(echo "@")
header_array=("${header_array[@]/$delete}")

##finding array keys 
for i in "${!header_array[@]}"; do 
  printf "%s\t%s\n" "$i" "${header_array[$i]}"
done

#begin assigning variables, follow this format
general_Fastq_file_header="$(printf '%s' '@instrument_id:run_number:flowcell_ID:lane:tile:x_pos:y_pos read:is_filtered:control_number:index_sequence')"
#FLOWCELL_BARCODE = @(instrument id):(run number):(flowcell ID)
IFS=': ' read -a general_Fastq_file_header_array <<< $general_Fastq_file_header
printf "%s\n" ${general_Fastq_file_header_array[@]}

#see fastq details
for i in "${!header_array[@]}"; do 
  printf "%s\t%s\t%s\n" "$i" "${header_array[$i]}" "${general_Fastq_file_header_array[$i]}" 
done
 
 #https://angus.readthedocs.io/en/2017/Read_group_info.html
# https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups


file_Path=$file
file_Name_on_drive=$(basename -s .fastq.gz $file)
SM=$(basename -s .fastq.gz $file | cut -d '_' -f 1) ##please remember to trim away any "R1/R2" to match sample name on terra tables, I assume file names are sperated by delimiter '_'

lane=${header_array[3]}
LB=$SM.${header_array[10]}
read=${header_array[7]}

#Sample_barcode=$(cat "tmp_header.txt" | head -n 1 | cut -d " " -f 2 | cut -d ":" -f 4)
flowcell_barcode=${header_array[0]}.${header_array[1]}.${header_array[2]}

PU=$flowcell_barcode.${header_array[3]}.${header_array[10]}
RGID=$flowcell_barcode.${header_array[3]}
#dump_RG_information_to_file=$(echo -e "$SM" ' \t ' "$read" ' \t ' "$PU" ' \t ' $RGID ' \t ' $LB ' \t ' $lane ' \t ' $file_Path ' \t ' ${header_array[10]})
echo -e "$SM" ' \t ' "$read" ' \t ' "$PU" ' \t ' "$RGID" ' \t ' "$LB" ' \t ' "$lane" ' \t ' "$file_Path" ' \t ' "${header_array[10]}" >> ${read_groups_all}.txt
done


