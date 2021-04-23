# Author; Samuel Ahuno
# Date; February, 25th 2021
# Purpose: cfDNA Anal cancer analysis- Get paramters

library(reshape2)
library(hablar)
library(lubridate)
library(tidyverse)
library(stringr)
library(data.table)




##function to extract ichorCNA parameters

IchorCNA_param_extract <- function(param_text_file=""){
  ###read param files
  read.ichor.param <- fread(param_text_file, header=FALSE,  col.names = c("parameter","value"),
  sep = "\t", strip.white = TRUE,blank.lines.skip = TRUE, nrows = 9,skip = 5)
  
  ###get basenames of files 
  sampleName <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(param_text_file)))
  
  ##remove any non-unifrm characters
  read.ichor.param$sample <- sampleName
  read.ichor.param$parameter <- gsub(" ","_",gsub(":","",read.ichor.param$parameter))
  
  ###make data wide
  ichor.params <- reshape2::dcast(read.ichor.param, sample ~ parameter, value.var = "value")
  ichor.params.final <- ichor.params %>% dplyr::select("sample","Tumor_Fraction","Ploidy",
  "Subclone_Fraction","Fraction_CNA_Subclonal","Fraction_Genome_Subclonal","Gender",
  "ChrX_median_log_ratio","ChrY_coverage_fraction","Coverage")

  return(ichor.params.final)
}


##create empty vector to hold results from loop
ichor.param.dataStore <- data.frame(stringsAsFactors=FALSE)

### for loop to run all samples
for (params_txt in list.files(path="/path/to/ichorCNA_params_text_files/",full.names=TRUE,pattern = "params.txt", recursive = TRUE)) {
  df.data <- IchorCNA_param_extract(params_txt)
  ichor.param.dataStore <- as.data.frame(rbind(ichor.param.dataStore,df.data))
}


#save a copy, just in case you miss a step
ichor.param.dataStore.copy <- ichor.param.dataStore

##coerce columns to numeric 
ichor.param.dataStore.copy <- ichor.param.dataStore.copy %>% 
hablar::convert(num(Tumor_Fraction, Ploidy,Subclone_Fraction,Fraction_CNA_Subclonal,
Fraction_Genome_Subclonal,ChrX_median_log_ratio,
ChrY_coverage_fraction,Coverage),chr(Gender))


##function to convert to perentage nicely
cal_percent <- function(x, na.rm=FALSE) (x * 100)
ichor.param.dataStore.copy <- ichor.param.dataStore.copy %>% 
mutate_at(vars(c("Tumor_Fraction","Subclone_Fraction","Fraction_CNA_Subclonal",
"Fraction_Genome_Subclonal","ChrY_coverage_fraction")),cal_percent)


#rename samples, to reflect genomic ids, files sometimes come with long basenames
# ichor.param.dataStore.copy$sample <- gsub("RP-2258_","",gsub("_v1_WGS_OnPrem","",ichor.param.dataStore.copy$sample))
# dim(ichor.param.dataStore.copy)

#optional write to files
write.table(ichor.param.dataStore.copy,file = "ichor_params_combined.tsv",quote=FALSE, sep='\t',row.names = F)