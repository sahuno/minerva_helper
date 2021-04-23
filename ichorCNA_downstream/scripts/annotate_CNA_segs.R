# Author; Samuel Ahuno
# Date; February, 25th 2021
# Purpose: cfDNA Anal cancer analysis

library(biomaRt)
library(GenomeInfoDb)
library(GenomicFeatures)
library(data.table)
library(data.table)
library(tidyverse)
library(reshape2)
library(dplyr)

##set path to script that annotates Copy number regions with human ref genes 
path_to_funAnnotateR <- "/Users/samuelahuno/Polak_lab10082019/Collabs_internal/analPreCancer_cfDNA_02252021/scripts/functionAnnotate.R"
annotated.seg <- paste0("../data/output/melt.df.annotations.genes.anal.tsv")

if (!file.exists(annotated.seg)) {
  #Path of ichorCNA seg files
  cna.seg.paths <- "/Users/samuelahuno/Polak_lab10082019/Collabs_internal/analPreCancer_cfDNA_02252021/data/raw_data/cna_seg"
  source(path_to_funAnnotateR ) ###load script to extract copy number from *cna.seg
  
  dataAnal <- data.frame(stringsAsFactors=FALSE)
  dat = list.files(path=cna.seg.paths,full.names=TRUE,pattern = "cna.seg", recursive = TRUE)[1]
  ### use this to annotate all samples together
  for (dat in list.files(path=cna.seg.paths,full.names=TRUE,pattern = "cna.seg", recursive = TRUE)) {
    df.data <- AnnotateEnsembleGenes.anal(geneList=gr.ensemblGenes.biomart.hg19,filePath=dat)
    dataAnal <- as.data.frame(c(dataAnal,df.data))
  }
  
  ##get gene names
  getGeneName<- dataAnal[,grepl("Gene",names(dataAnal))][1]
  getCytobands <- dataAnal[,grepl("cytoband",names(dataAnal))][1]
  
  removeMatch <- c("Gene","cytoband") ### patterns to remove
  ##get cna nunbers
  df.annotations.genes <- cbind(getGeneName,getCytobands,dataAnal[,!grepl(paste(removeMatch, collapse = "|"),names(dataAnal))])
  
  ###insert NA's into empty cells and get rid of all empty cases get
  df.annotations.genes[df.annotations.genes==""]<-NA
  df.annotations.genes <- df.annotations.genes[complete.cases(df.annotations.genes), ]
  
  ##get unique  cells 
  df.annotations.genes <- distinct(df.annotations.genes, Gene_Names, .keep_all = TRUE)
  
  
  ##add band to gene names
  df.annotations.genes.loci  <- df.annotations.genes 
  df.annotations.genes.loci$gene_loci <- NA
  df.annotations.genes.loci$geneBand <- paste0(df.annotations.genes.loci$Gene_Names," (",df.annotations.genes.loci$cytoband,")")
  df.saved.gene.loci <- df.annotations.genes.loci[,c(1,2,ncol(df.annotations.genes.loci))] ### save gene name and locis, use later to interect
  
  ##sanity check #ping!
  df.annotations.genes %>% filter (Gene_Names == "ERBB2")
  df.annotations.genes %>% filter (Gene_Names == "CD274")
  
  rownames(df.annotations.genes) <- df.annotations.genes[,1] 
  df.annotations.genes.final <- df.annotations.genes[-1]
  
  ##comut ready data - for plotting
  melt.df.annotations.genes.anal <-  melt(df.annotations.genes[-2],id="Gene_Names",id.name="Gene",variable.name = "sample",value.name ="corrected_copy_number")
  
  ##write annotated segments to file
  write.table(melt.df.annotations.genes.anal,file = "../data/output/melt.df.annotations.genes.anal.tsv",quote=FALSE, sep='\t',row.names = F)
  } else{
  melt.df.annotations.genes.anal <- read.delim(annotated.seg, sep = "\t", header = TRUE, quote = "\"", dec = ".", fill = TRUE, comment.char = "#") ##read files
}







##subset selected genes for comut
anal_cancer_specific_genes <- str_sort(c("CCL22","STK11","DNMT3B","TERC","PIK3CA","TP63","TGFBR2","TERT","DUSP4","ATM","BAP1","TRAF3","TSC1","CD274"), numeric = TRUE)
sub.melt.df.annotations.cna.genes.anal = melt.df.annotations.genes.anal[which(melt.df.annotations.genes.anal$Gene_Names %in% anal_cancer_specific_genes),]

melt.df.annotations.genes.anal %>% filter(corrected_copy_number < 1)

###sumamry stats on 
sub.melt.df.annotations.cna.genes.summary.anal <- sub.melt.df.annotations.cna.genes.anal %>% group_by(Gene_Names) %>% 
  dplyr::summarise(count=n(),minCN=min(corrected_copy_number),
                                                              medianCN=median(corrected_copy_number),maxCN=max(corrected_copy_number),
                                                              CNgain=sum(corrected_copy_number>2 & corrected_copy_number<4),
                                                              CNgain_Percent=((sum(corrected_copy_number>2 & corrected_copy_number<4))/n())*100,
                                                              CNamp_Percent=((sum(corrected_copy_number==4))/n())*100,
                                                              CNhamp_Percent=((sum(corrected_copy_number>4))/n())*100,
                                                              CNneutr_Percent=((sum(corrected_copy_number==2))/n())*100,
                                                              CNdel_Perc=((sum(corrected_copy_number<2))/n())*100,
                                                              CN_anyCN_above1_Percent=((sum(corrected_copy_number>1))/n())*100,
                                                              CN_atLeast_gain=(sum(corrected_copy_number>2)),
                                                               CN_atLeast_gain_Percent=((sum(corrected_copy_number>2))/n())*100)

write.table(sub.melt.df.annotations.cna.genes.summary.anal, file="../data/output/statistics_copyNumber.tsv",quote=FALSE,row.names = FALSE,sep='\t')
###plot percentage of individuals with any cn higher tan neutral in the cohort
#ggplot(sub.melt.df.annotations.cna.genes.summary.anal, aes(x=as.factor(Gene_Names),y=CN_anyCN_above1_Percent)) + geom_col(color="blue")



sub.melt.df.annotations.cna.genes.anal$cnaEvent <- NA
test.cna.event.anal <- sub.melt.df.annotations.cna.genes.anal %>% mutate(cnaEvent = ifelse(corrected_copy_number==2,"Neutral",ifelse(corrected_copy_number == 1, "Homodeletion",ifelse(corrected_copy_number == 4, "Amplification",ifelse(corrected_copy_number > 4, "High amplification",
                                                                                                                                                                                                                                          ifelse(corrected_copy_number == 3, "Gain", cnaEvent))))))

##write to file
test.cna.event.anal.cmut.renamed <- test.cna.event.anal
test.cna.event.anal.cmut.renamed$sample <- gsub("RP.2258_","",gsub("_v1_WGS_OnPrem","",test.cna.event.anal.cmut.renamed$sample))
test.cna.event.anal.cmut.renamed <- test.cna.event.anal.cmut.renamed[,c(2,1,4)]
  
names(test.cna.event.anal.cmut.renamed) <- c("sample", "category","value")
#write.table(test.cna.event.anal.cmut.renamed, file="../data/output/coMutPlot_ready_data_allSamples.anal.tsv",quote=FALSE, col.names = c("sample", "category","value"),sep='\t',row.names = FALSE)


##read pilot data
path.pilot.cna <- "/Users/samuelahuno/Polak_lab10082019/Collabs_internal/analPreCancer_ULPS_pilot/data/output/coMutPlot_ready_data_allSamples.anal.tsv"
cna.anal.pilot <- read.delim(path.pilot.cna, sep = "\t", header = TRUE, quote = "\"", dec = ".", fill = TRUE, comment.char = "#") ##read files
hsil_ain3 <- c("HSIL_P06","HSIL_P01","HSIL_P04",  "HSIL_P07") #for pilot data #c("BRP23739","BRP23844","BRP23737","BRP23845")

df.cna.anal.pilot.ain3 <- cna.anal.pilot %>% filter(sample %in% c("BRP23739","BRP23844","BRP23737","BRP23845"))
df.old.new.comut <- rbind(test.cna.event.anal.cmut.renamed,df.cna.anal.pilot.ain3)

write.table(df.old.new.comut, file="../data/output/coMutPlot_ready_data_pilot_new_samples.anal.tsv",quote=FALSE, col.names = c("sample", "category","value"),sep='\t',row.names = FALSE)


