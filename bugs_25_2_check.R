## PUBRE 
##############
rm(list = ls())
library(dplyr)
library(tidyr)

library(data.table)
library(arrow)

papers <- fread("/home/janet/Documents/disgenetplus_v2/DATA/papers_disgenet_v25_2_bugs.tsv")
papers

# data <- fread("/home/janet/Downloads/prod_dgn_pubre_v25.3-dev-tmvar3_part_1.tsv")
data <- fread("/home/janet/Downloads/prod_dgn_pubre_v25.3-dev-tmvar2_part_1.tsv")
data <- data[ id %in% papers$x, ]

length(intersect(papers$x, data$id))
length(setdiff(papers$x, data$id))
# 148 154
data

pubre <- data %>% unique()

geneinfo <- fread("~/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/CURATED_SOURCES/CURATED_SOURCES_DOWNLOADS/Homo_sapiens.gene_info.gz")
geneinfo <- geneinfo %>% select(GeneID, Symbol)
genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/SHARED/PERCOLATORS/v25.1.1/GENE/TSV/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv" )
disnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/original_data/UMLS/dis_names_2024AB.tsv", quote = "") 
chemnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/CHEMICAL_PERCOLATOR/chemical_percolator_names_25.3.tsv", sep = "\t",quote="")


pubre$global_gene_ID <- gsub("\\'", "", pubre$e1_gene_ids)
pubre$global_gene_ID <- gsub("\\[", "", pubre$global_gene_ID)
pubre$global_gene_ID <- gsub("\\]", "", pubre$global_gene_ID)
pubre$global_gene_ID <- gsub("-", "_", pubre$global_gene_ID)
setdiff(pubre$global_gene_ID, genes$global_gene_ID)

pubre <- separate_rows(pubre, global_gene_ID , sep = ", ")
setdiff(pubre$global_gene_ID, genes$global_gene_ID)


pubre$cui <- gsub("\\'", "", pubre$e2_dis_cuis)
pubre$cui <- gsub("\\[", "", pubre$cui)
pubre$cui <- gsub("\\]", "", pubre$cui)
pubre <- separate_rows(pubre, cui , sep = ", ")
setdiff(pubre$cui, disnames$V1)

pubre <- merge(pubre, genes,by = "global_gene_ID", all.x = T)
pubre <- pubre %>% select(-UNIPROT, -HGNC)
setdiff(pubre$NCBI, geneinfo$GeneID)
pubre <- separate_rows(pubre, NCBI , sep = ", ")
setdiff(pubre$NCBI, geneinfo$GeneID)

pubre <- merge(pubre, geneinfo,by.x = "NCBI", by.y = "GeneID", all.x=T )

pubre$chemical_ids  <- pubre$e3_chemical_ids
pubre$chemical_ids <- gsub("\\'", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("\\[", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("\\]", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("-", "_", pubre$chemical_ids)

setdiff(pubre$chemical_ids, chemnames$concept_id)
pubre <- separate_rows(pubre, chemical_ids , sep = ", ")

setdiff(pubre$chemical_ids, chemnames$concept_id)

pubre[! pubre$chemical_ids %in% chemnames$concept_id ,] %>% group_by(e3_text,chemical_ids) %>% summarise(nn = n()) %>% arrange(desc(nn))

pubre <- merge(pubre, chemnames, by.x = "chemical_ids", by.y =  "concept_id", all.x = T)
pubre %>% group_by(e3_text,chemical_ids) %>% summarise(nn = n()) %>% arrange(desc(nn))
pubre <- merge(pubre, disnames, by.x = "cui", by.y =  "V1", all.x = T)

write.table(pubre, "/home/janet/Documents/disgenetplus_v2/DATA/cleaned_pubre_v25_3_tmvar2.tsv",
            quote = F, row.names= F, sep = "\t" )










## PUBRE 25.2
##############
rm(list = ls())
library(dplyr)
library(tidyr)

library(data.table)
library(arrow)

papers <- fread("/home/janet/Documents/disgenetplus_v2/DATA/papers_disgenet_v25_2_bugs.tsv")
papers

dir <- "/media/janet/Janet/disgenetplus_v2/DATA/pubre_v25_2"
dir <- "~/Documents/disgenetplus_v2/DATA/pubre_v25_2"
files <- list.files(dir)
files <- files[ grep("prod_dgn_pubre_v25.", files)]
filename <- paste(dir,  files[1], sep = "/")
data <- fread(filename)
data <- data[ id %in% papers$x, ]

for (ff in files[2:length(files)]){
  filename <- paste(dir, ff, sep = "/")
  print(filename)
  if (file.exists(filename)){
    tmp <- fread(filename, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE   )
    tmp <- tmp[ id %in% papers$x,]
    if(nrow(tmp )> 0){
    data<- rbind(data, tmp)
    }
      data<- rbind(data, tmp)
  }
  else {
    print("file does not exit")
  }
}


length(intersect(papers$x, data$id))
length(setdiff(papers$x, data$id))


data
pubre <- data %>% unique()
# rm(data)


geneinfo <- fread("~/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/CURATED_SOURCES/CURATED_SOURCES_DOWNLOADS/Homo_sapiens.gene_info.gz")
geneinfo <- geneinfo %>% select(GeneID, Symbol)
genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/SHARED/PERCOLATORS/v25.1.1/GENE/TSV/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv" )
#genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/NLP/NER_CREATION/DATA/nerGeneCreate/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv")
disnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/original_data/UMLS/dis_names_2024AB.tsv", quote = "") 
chemnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/CHEMICAL_PERCOLATOR/chemical_percolator_names_25.1.tsv", sep = "\t",quote="")


pubre$global_gene_ID <- gsub("\\'", "", pubre$e1_gene_ids)
pubre$global_gene_ID <- gsub("\\[", "", pubre$global_gene_ID)
pubre$global_gene_ID <- gsub("\\]", "", pubre$global_gene_ID)
pubre$global_gene_ID <- gsub("-", "_", pubre$global_gene_ID)
setdiff(pubre$global_gene_ID, genes$global_gene_ID)

pubre <- separate_rows(pubre, global_gene_ID , sep = ", ")
setdiff(pubre$global_gene_ID, genes$global_gene_ID)


pubre$cui <- gsub("\\'", "", pubre$e2_dis_cuis)
pubre$cui <- gsub("\\[", "", pubre$cui)
pubre$cui <- gsub("\\]", "", pubre$cui)
pubre <- separate_rows(pubre, cui , sep = ", ")
setdiff(pubre$cui, disnames$V1)

#pubre <- merge(pubre, disnames, by.x = "cui", by.y =  "V1", all.x = T)
# pubre <- merge(pubre, disnames, by.x = "cui", by.y =  "V1" )

pubre <- merge(pubre, genes,by = "global_gene_ID", all.x = T)
pubre <- pubre %>% select(-UNIPROT, -HGNC)
setdiff(pubre$NCBI, geneinfo$GeneID)
pubre <- separate_rows(pubre, NCBI , sep = ", ")
setdiff(pubre$NCBI, geneinfo$GeneID)

# pubre <- merge(pubre, geneinfo,by.x = "NCBI", by.y = "GeneID", all.x=T )

pubre$chemical_ids  <- pubre$e3_chemical_ids
pubre$chemical_ids <- gsub("\\'", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("\\[", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("\\]", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("-", "_", pubre$chemical_ids)

setdiff(pubre$chemical_ids, chemnames$concept_id)
pubre <- separate_rows(pubre, chemical_ids , sep = ", ")

setdiff(pubre$chemical_ids, chemnames$concept_id)

pubre[! pubre$chemical_ids %in% chemnames$concept_id ,] %>% group_by(e3_text,chemical_ids) %>% summarise(nn = n()) %>% arrange(desc(nn))

pubre <- merge(pubre, chemnames, by.x = "chemical_ids", by.y =  "concept_id", all.x = T)
pubre %>% group_by(e3_text,chemical_ids) %>% summarise(nn = n()) %>% arrange(desc(nn))


pubre <- merge(pubre, geneinfo,by.x = "NCBI", by.y = "GeneID", all.x = T)
setdiff(pubre$cui, disnames$V1)

pubre <- merge(pubre, disnames, by.x = "cui", by.y =  "V1", all.x = T)
# pubre <- merge(pubre, disnames, by.x = "cui", by.y =  "V1" )


write.table(pubre, "/home/janet/Documents/disgenetplus_v2/DATA/cleaned_pubre_v25_2_bug.tsv",
            quote = F, row.names= F, sep = "\t" )


## PUBNER 25.2
##############
rm(list = ls())
library(dplyr)
library(tidyr)

library(data.table)
library(arrow)



geneinfo <- fread("~/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/CURATED_SOURCES/CURATED_SOURCES_DOWNLOADS/Homo_sapiens.gene_info.gz")
geneinfo <- geneinfo %>% select(GeneID, Symbol)
genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/SHARED/PERCOLATORS/v25.1.1/GENE/TSV/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv" )
#genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/NLP/NER_CREATION/DATA/nerGeneCreate/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv")
disnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/original_data/UMLS/dis_names_2024AB.tsv", quote = "") 
chemnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/CHEMICAL_PERCOLATOR/chemical_percolator_names_25.1.tsv", sep = "\t",quote="")


papers <- fread("/home/janet/Documents/disgenetplus_v2/DATA/papers_disgenet_v25_2_bugs.tsv")
papers

dir <- "/media/janet/Janet/disgenetplus_v2/DATA/pubner_v25_2/"
dir <- "~/Documents/disgenetplus_v2/DATA/pubner_v25_2"
files <- list.files(dir)
files <- files[ grep("prod_dgn_pubner_v25.", files)]

filename <- paste(dir,  paste0("prod_dgn_pubner_v25.2_part_1.tsv"), sep = "/")
data <- fread(filename, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE   )
data <- data[ id %in% papers$x,]

for (ff in 2:length(files)){
  filename <- paste(dir,  paste0("prod_dgn_pubner_v25.2_part_", ff,".tsv"), sep = "/")
  # print(filename)
  if (file.exists(filename)){
    
    tmp <- fread(filename, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE   )
    tmp<- tmp[id %in% papers$x,]
    data<- rbind(data, tmp)
  }
  else {
    print("file does not exit")
  }
}

pubner <- data
length(setdiff(papers$x, pubner$id))


pubner$global_gene_ID <- gsub("\\'", "", pubner$gene_ids)
pubner$global_gene_ID <- gsub("\\[", "", pubner$global_gene_ID)
pubner$global_gene_ID <- gsub("\\]", "", pubner$global_gene_ID)
pubner$global_gene_ID <- gsub("-", "_", pubner$global_gene_ID)
setdiff(pubner$global_gene_ID, genes$global_gene_ID)

pubner <- separate_rows(pubner, global_gene_ID , sep = ", ")
setdiff(pubner$global_gene_ID, genes$global_gene_ID)

pubner$cui <- gsub("\\'", "", pubner$dis_cuis)
pubner$cui <- gsub("\\[", "", pubner$cui)
pubner$cui <- gsub("\\]", "", pubner$cui)

setdiff(pubner$cui, disnames$V1)
pubner <- separate_rows(pubner, cui , sep = ", ")
setdiff(pubner$cui, disnames$V1)

pubner <- merge(pubner, disnames, by.x = "cui", by.y =  "V1", all.x = T)

# pubner <- pubner[ ! (pubner$type == "variant" & pubner$var_snp_list == "[]" ), ]
# pubner <- pubner[ pubner$type != "negation_marker",]
pubner <- merge(pubner, genes,by = "global_gene_ID", all.x = T)
pubner <- pubner %>% select(-UNIPROT, -HGNC)

# pubner <- pubner %>% select(id, section, type, text, start, end, NCBI, var_snp_list, cui, V15, faiss_norm_score,ner_model_score) %>% unique()
pubner <- merge(pubner, geneinfo,by.x = "NCBI", by.y = "GeneID", all.x = T)

pubner$original_chemical_ids  <- pubner$chemical_ids 
pubner$chemical_ids <- gsub("\\'", "", pubner$chemical_ids)
pubner$chemical_ids <- gsub("\\[", "", pubner$chemical_ids)
pubner$chemical_ids <- gsub("\\]", "", pubner$chemical_ids)
pubner$chemical_ids <- gsub("-", "_", pubner$chemical_ids)

setdiff(pubner$chemical_ids, chemnames$concept_id)
pubner <- separate_rows(pubner, chemical_ids , sep = ", ")

setdiff(pubner$chemical_ids, chemnames$concept_id)

pubner[! pubner$chemical_ids %in% chemnames$concept_id & pubner$type == "chemical",] %>% group_by(text,chemical_ids) %>% 
  summarise(nn = n()) %>% arrange(desc(nn)) %>% print(n = 200)

pubner <- merge(pubner, chemnames, by.x = "chemical_ids", by.y =  "concept_id", all.x = T)
pubner[ pubner$chemical_ids %in% chemnames$concept_id & pubner$type == "chemical",] %>% group_by(text,chemical_ids, chemical_name) %>% 
  summarise(nn = n()) %>% arrange(desc(nn)) %>% print(n = 200)

# pubner <- pubner %>% rename(pmid = id, dbsnpid = var_snp_list, diseaseid = cui, disease_name = V15, gene_symbol = Symbol)
# pubner <- pubner %>% select(NCBI, gene_symbol,  type, text, start, end, diseaseid, disease_name, pmid, section )


write.table(pubner, "/home/janet/Documents/disgenetplus_v2/DATA/cleaned_pubner_v25_2_bug.tsv",
            quote = F, row.names= F, sep = "\t" )

## PUBRE 25.3
##############
rm(list = ls())
library(dplyr)
library(tidyr)

library(data.table)
library(arrow)

papers <- fread("/home/janet/Documents/disgenetplus_v2/DATA/papers_disgenet_v25_2_bugs.tsv")
papers

dir <- "/media/janet/Janet/disgenetplus_v2/DATA/pubre_v25_3_dev/"
dir <- "/home/janet/Downloads/"
files <- list.files(dir)
files <- files[ grep("prod_dgn_pubre_v25.", files)]
filename <- paste(dir,  files[1], sep = "/")
data <- fread(filename)
data <- data[ id %in% papers$x, ]

for (ff in files[2:length(files)]){
  filename <- paste(dir, ff, sep = "/")
  print(filename)
  if (file.exists(filename)){
    tmp <- fread(filename, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE   )
    tmp <- tmp[ id %in% papers$x,]
    if(nrow(tmp )> 0){
      data<- rbind(data, tmp)
    }
    data<- rbind(data, tmp)
  }
  else {
    print("file does not exit")
  }
}


length(intersect(papers$x, data$id))
length(setdiff(papers$x, data$id))


data

pubre <- data %>% unique()
# rm(data)


geneinfo <- fread("~/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/CURATED_SOURCES/CURATED_SOURCES_DOWNLOADS/Homo_sapiens.gene_info.gz")
geneinfo <- geneinfo %>% select(GeneID, Symbol)
genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/SHARED/PERCOLATORS/v25.1.1/GENE/TSV/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv" )
#genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/NLP/NER_CREATION/DATA/nerGeneCreate/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv")
disnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/original_data/UMLS/dis_names_2024AB.tsv", quote = "") 
chemnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/CHEMICAL_PERCOLATOR/chemical_percolator_names_25.3.tsv", sep = "\t",quote="")


pubre$global_gene_ID <- gsub("\\'", "", pubre$e1_gene_ids)
pubre$global_gene_ID <- gsub("\\[", "", pubre$global_gene_ID)
pubre$global_gene_ID <- gsub("\\]", "", pubre$global_gene_ID)
pubre$global_gene_ID <- gsub("-", "_", pubre$global_gene_ID)
setdiff(pubre$global_gene_ID, genes$global_gene_ID)

pubre <- separate_rows(pubre, global_gene_ID , sep = ", ")
setdiff(pubre$global_gene_ID, genes$global_gene_ID)


pubre$cui <- gsub("\\'", "", pubre$e2_dis_cuis)
pubre$cui <- gsub("\\[", "", pubre$cui)
pubre$cui <- gsub("\\]", "", pubre$cui)
pubre <- separate_rows(pubre, cui , sep = ", ")
setdiff(pubre$cui, disnames$V1)

pubre <- merge(pubre, genes,by = "global_gene_ID", all.x = T)
pubre <- pubre %>% select(-UNIPROT, -HGNC)
setdiff(pubre$NCBI, geneinfo$GeneID)
pubre <- separate_rows(pubre, NCBI , sep = ", ")
setdiff(pubre$NCBI, geneinfo$GeneID)

 pubre <- merge(pubre, geneinfo,by.x = "NCBI", by.y = "GeneID", all.x=T )

pubre$chemical_ids  <- pubre$e3_chemical_ids
pubre$chemical_ids <- gsub("\\'", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("\\[", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("\\]", "", pubre$chemical_ids)
pubre$chemical_ids <- gsub("-", "_", pubre$chemical_ids)

setdiff(pubre$chemical_ids, chemnames$concept_id)
pubre <- separate_rows(pubre, chemical_ids , sep = ", ")

setdiff(pubre$chemical_ids, chemnames$concept_id)

pubre[! pubre$chemical_ids %in% chemnames$concept_id ,] %>% group_by(e3_text,chemical_ids) %>% summarise(nn = n()) %>% arrange(desc(nn))

pubre <- merge(pubre, chemnames, by.x = "chemical_ids", by.y =  "concept_id", all.x = T)
pubre %>% group_by(e3_text,chemical_ids) %>% summarise(nn = n()) %>% arrange(desc(nn))
pubre <- merge(pubre, disnames, by.x = "cui", by.y =  "V1", all.x = T)

write.table(pubre, "/home/janet/Documents/disgenetplus_v2/DATA/cleaned_pubre_v25_3_dev_bug.tsv",
            quote = F, row.names= F, sep = "\t" )



## PUBNER 25.3
##############
rm(list = ls())
library(dplyr)
library(tidyr)

library(data.table)
library(arrow)



geneinfo <- fread("~/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/CURATED_SOURCES/CURATED_SOURCES_DOWNLOADS/Homo_sapiens.gene_info.gz")
geneinfo <- geneinfo %>% select(GeneID, Symbol)
genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/SHARED/PERCOLATORS/v25.1.1/GENE/TSV/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv" )
#genes <- fread("/home/janet/Documents/disgenetplus_v2/DISGENET_PLUS_RESOURCES/NLP/NER_CREATION/DATA/nerGeneCreate/disgenet_gene_dictionary_global_to_specific_gene_ID_mappings_term_added.tsv")
disnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/original_data/UMLS/dis_names_2024AB.tsv", quote = "") 
chemnames <- fread("/home/janet/Documents/disgenetplus_v2/DATA/CHEMICAL_PERCOLATOR/chemical_percolator_names_25.3.tsv", sep = "\t",quote="")


papers <- fread("/home/janet/Documents/disgenetplus_v2/DATA/papers_disgenet_v25_2_bugs.tsv")
papers

dir <- "/media/janet/Janet/disgenetplus_v2/DATA/pubner_v25_3_dev/"
dir <- "/home/janet/Downloads/"
files <- list.files(dir)
files <- files[ grep("prod_dgn_pubner_v25.", files)]

filename <- paste(dir,  files[1], sep = "/")
data <- fread(filename, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE   )
data <- data[ id %in% papers$x,]

for (ff in 2:length(files)){
  filename <- paste(dir,  paste0("prod_dgn_pubner_v25.3-dev_part_", ff,".tsv"), sep = "/")
  # print(filename)
  if (file.exists(filename)){
    
    tmp <- fread(filename, sep = "\t", header = T, stringsAsFactors = F, fill=TRUE   )
    tmp<- tmp[id %in% papers$x,]
    data<- rbind(data, tmp)
  }
  else {
    print("file does not exit")
  }
}

pubner <- data
length(setdiff(papers$x, pubner$id))


pubner$global_gene_ID <- gsub("\\'", "", pubner$gene_ids)
pubner$global_gene_ID <- gsub("\\[", "", pubner$global_gene_ID)
pubner$global_gene_ID <- gsub("\\]", "", pubner$global_gene_ID)
pubner$global_gene_ID <- gsub("-", "_", pubner$global_gene_ID)
setdiff(pubner$global_gene_ID, genes$global_gene_ID)

pubner <- separate_rows(pubner, global_gene_ID , sep = ", ")
setdiff(pubner$global_gene_ID, genes$global_gene_ID)

pubner$cui <- gsub("\\'", "", pubner$dis_cuis)
pubner$cui <- gsub("\\[", "", pubner$cui)
pubner$cui <- gsub("\\]", "", pubner$cui)

setdiff(pubner$cui, disnames$V1)
pubner <- separate_rows(pubner, cui , sep = ", ")
setdiff(pubner$cui, disnames$V1)

pubner <- merge(pubner, disnames, by.x = "cui", by.y =  "V1", all.x = T)

# pubner <- pubner[ ! (pubner$type == "variant" & pubner$var_snp_list == "[]" ), ]
# pubner <- pubner[ pubner$type != "negation_marker",]
pubner <- merge(pubner, genes,by = "global_gene_ID", all.x = T)
pubner <- pubner %>% select(-UNIPROT, -HGNC)

# pubner <- pubner %>% select(id, section, type, text, start, end, NCBI, var_snp_list, cui, V15, faiss_norm_score,ner_model_score) %>% unique()
pubner <- merge(pubner, geneinfo,by.x = "NCBI", by.y = "GeneID", all.x = T)

pubner$chemical_ids_original <- pubner$chemical_ids 
pubner$chemical_ids <- gsub("\\'", "", pubner$chemical_ids)
pubner$chemical_ids <- gsub("\\[", "", pubner$chemical_ids)
pubner$chemical_ids <- gsub("\\]", "", pubner$chemical_ids)
pubner$chemical_ids <- gsub("-", "_", pubner$chemical_ids)

setdiff(pubner$chemical_ids, chemnames$concept_id)
pubner <- separate_rows(pubner, chemical_ids , sep = ", ")

setdiff(pubner$chemical_ids, chemnames$concept_id)

pubner[! pubner$chemical_ids %in% chemnames$concept_id & pubner$type == "chemical",] %>% group_by(text,chemical_ids) %>% 
  summarise(nn = n()) %>% arrange(desc(nn)) %>% print(n = 200)

pubner <- merge(pubner, chemnames, by.x = "chemical_ids", by.y =  "concept_id", all.x = T)

pubner[ pubner$chemical_ids %in% chemnames$concept_id & pubner$type == "chemical",] %>% group_by(text,chemical_ids, chemical_name) %>% 
  summarise(nn = n()) %>% arrange(desc(nn)) %>% print(n = 200)



# pubner <- pubner %>% rename(pmid = id, dbsnpid = var_snp_list, diseaseid = cui, disease_name = V15, gene_symbol = Symbol)
# pubner <- pubner %>% select(NCBI, gene_symbol,  type, text, start, end, diseaseid, disease_name, pmid, section )


write.table(pubner, "/home/janet/Documents/disgenetplus_v2/DATA/cleaned_pubner_v25_3_bug.tsv",
            quote = F, row.names= F, sep = "\t" )
