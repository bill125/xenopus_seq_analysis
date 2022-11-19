library(WGCNA)
library(reshape2)
library(tidyverse)

#Note: This example performs QENIE on the liver-to-adipose circuit.
#Preproccesing of data from the original study (GSE64770) is listed in the README.md file
gene_analysis <- function(x, y) {
  # Import secreted peptides
  Secreted_proteins <- read.delim("data/xen_secreted_protein.csv", sep=",", header=T, check.names=F)
  
  # Import your data
  x_data <- read.delim(paste("data/input/", x, "_vs_", y, "-", x, ".csv", sep=""), check.names=F)
  y_data <- read.delim(paste("data/input/", x, "_vs_", y, "-", y, ".csv", sep=""), check.names=F)
  # Using these two datasets we will proceed:
  
  # Construct cross-tissue correlation and pvalue matrices
  tmp = bicorAndPvalue(x_data, y_data, use='pairwise.complete.obs')
  
  # Compute rowsum -log(pvalue) for ranking interactions and filter for factors proteins as secreted
  scores = rowSums(-log(tmp$p), na.rm=TRUE)
  
  # Retain only secreted peptides, normalize the Ssec score by number of target tissue probes (In this case, adipose contains 12242 genes)
  # order by "sig score"
  Ssec = data.frame(Gene_symbol = names(scores), score = scores) %>%
    filter(Gene_symbol %in% Secreted_proteins$`Gene_names`) %>%
    mutate(Ssec = score / length(colnames(y_data))) %>%
    select(-score) %>%
    arrange(desc(Ssec))
  
  # return(score)
  
  write.table(Ssec, file=paste("result/", x, " X ", y, " ranked by sig score", sep=""), row.names=F, col.names=T, sep='\t', quote=F)
  
  #This produces a table to each secreted protein and its respective significance score across adipose transcripts
  #Note that Notum is listed as the 5th with an Sssec of 4.148
  #For our pipeline, we next check the tissue-specificty using BioGPS - note that this step is not necessary, but makes us more confident when conditioning the pathway enrichment.  This moves Notum up to the 2nd ranked protein
  
  #Condition correlation matrix on a by-gene basis for pathway enrichment - this example will focus on the protein, Notum
  #remove gene of interest (Notum) from the correlation matrix
  #the rownames of the file correspond to gene symbols for pathway enrichment, whereas the second column contains the bicor coefficent

  bicor.data = melt(tmp$bicor)
  colnames(bicor.data) = c('Gene_symbol_1', 'Gene_symbol_2', 'bicor')
  Notum = filter(bicor.data, Gene_symbol_1 == 'chga') %>%
    select(-Gene_symbol_1) %>%
    rename(Gene_symbol = Gene_symbol_2)
  
  #return(list("x_data"=x_data,"y_data"=y_data, "tmp"=tmp, "scores"=scores, "bicordata"=bicor.data))
  
  # Positively correlated pathways - "enhanced by protein"
  upreg = Notum %>% arrange(desc(bicor)) %>% head(500)
  write.table(upreg, file=paste("result/Positive Notum ", x, " X ", y, " Pathways Enrichment File", sep=""), col.names=F, sep='\t', quote=F)
  
  # Negatively correlated pathways - "suppressed by protein"
  downreg = Notum %>% arrange(bicor) %>% head(500)
  write.table(downreg, file=paste("result/Negative Notum ", x, " X ", y, " Pathways Enrichment File", sep=""), col.names=F, sep='\t', quote=F)
  
  # All pathways engaged by protein
  totalreg = Notum %>% arrange(desc(abs(bicor))) %>% head(500)
  write.table(totalreg, file=paste("result/Absolute Notum ", x, " X ", y, " Pathways Enrichment File", sep=""), col.names=F, sep='\t', quote=F)

}

organs = c("bra", "hea", "int", "liv", "lun", "mus", "ski")
ret = gene_analysis("int", "liv")
for (x in organs) {
  for (y in organs) {
    if (x != y) {
      break
      print(paste(x, y))
      gene_analysis(x, y)
    }   
  }
}
    
  
