#Comparison MA5
library(dplyr)

MA5_hmdb_input_MF<- read.xlsx("MA5Comparison/MA5InputMetaFetcheR.xlsx",sheetName = "discoveryHMDB")
MA5_hmdb_input_MF[c("kegg_id", "chebi_id","pubchem_id")] <- NA
MA5_kegg_input_MF<- read.xlsx("MA5Comparison/MA5InputMetaFetcheR.xlsx",sheetName = "discoveryKEGG")
MA5_kegg_input_MF[c("hmdb_id", "chebi_id","pubchem_id")] <- NA
MA5_pubchem_input_MF<- read.xlsx("MA5Comparison/MA5InputMetaFetcheR.xlsx",sheetName = "discoveryPubChem")
MA5_pubchem_input_MF[c("hmdb_id", "chebi_id","kegg_id")] <- NA

#--------------------Using Names to match MA5 and also for MetaFetcheR all IDS to map for Diamanti et al------------------------------------------
MA5_hmdb_output<- read.csv("MA5Comparison/MA_HMDB.csv")
MA5_kegg_output<- read.csv("MA5Comparison/MA_kegg.csv")
MA5_pubchem_output<- read.csv("MA5Comparison/MA_PubChem.csv")

# By names
query_results_tableMA5<-read.csv("MA5Comparison/MA_names.csv")
tempresultsKlevMA5<-resolve_metabolites(KlevIDSSelected)
tempresultsKlevMA5<-transformMFoutput(tempresultsKlevMA5$df,lipidmaps=FALSE)
rownames(query_results_tableMA5) <- query_results_tableMA5[,1]
query_results_tableMA5=as.data.frame(query_results_tableMA5)

NotFoundinMA5=which(!(trimws(KlevIds$common.name) %in% rownames(query_results_tableMA5)))
NotFoundinKlevMA5=which(!(rownames(query_results_tableMA5) %in% trimws(KlevIds$common.name)[-NotFoundinMA5]))
tempresultsKlevMA5=tempresultsKlevMA5[-NotFoundinMA5,]
tempKlevIDMA5=KlevIds[-NotFoundinMA5,]
tempresultsMA5=query_results_tableMA5[-NotFoundinKlevMA5,]
#order based on names in klevIDs
tempresultsMA5=tempresultsMA5[match(tempKlevIDMA5$common.name,rownames(tempresultsMA5)),]

tempresultsMA5=tempresultsMA5%>%
  mutate(across(everything(), as.character))
tempresultsKlevMA5=tempresultsKlevMA5 %>%
  mutate(across(everything(), as.character))
tempresultsMA5[is.na(tempresultsMA5)] <- "None"


#By IDS
MA5_kegg_output_cancer<- read.csv("MA5Comparison/MA_kegg_Priolo.csv")

#For priolo et al we can use the already processed data for comparison
#cancerDataMFresults

MA5_hmdb_output_MF<-resolve_metabolites(MA5_hmdb_input_MF)

#unlist everything in the output of MF

MA5_hmdb_output_MF<-transformMFoutput(MA5_hmdb_output_MF$df,lipidmaps=FALSE)

MA5_hmdb_output_MF=MA5_hmdb_output_MF %>%
  mutate(across(everything(), as.character))


MA5_kegg_output_MF<-resolve_metabolites(MA5_kegg_input_MF)

MA5_kegg_output_MF<-transformMFoutput(MA5_kegg_output_MF$df,lipidmaps=FALSE)

MA5_kegg_output_MF=MA5_kegg_output_MF %>%
  mutate(across(everything(), as.character))


MA5_pubchem_output_MF<-resolve_metabolites(MA5_pubchem_input_MF)

MA5_pubchem_output_MF<-transformMFoutput(MA5_pubchem_output_MF$df,lipidmaps =FALSE)

MA5_pubchem_output_MF=MA5_pubchem_output_MF %>%
  mutate(across(everything(), as.character))


query_results_tableMA5<-read.csv("MA5Comparison/MA_names.csv")

#---------MetaboAnalyst5.0and MetaFetcher using Diamanti et al data---------------------------------------

mapply(compareIDs,tempresultsMA$kegg_id,tempKlevID$KEGG)

length(which(mapply(compareIDs,tempresultsMA5$KEGG,tempresultsKlevMA5$kegg_id)))
length(which(mapply(compareIDs,tempresultsMA5$HMDB,tempresultsKlevMA5$hmdb_id)))
length(which(mapply(compareIDs,tempresultsMA5$ChEBI,tempresultsKlevMA5$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))
length(which(mapply(compareIDs,tempresultsMA5$PubChem,tempresultsKlevMA5$pubchem_id)))

#----MetaboAnalyst------------------
#NULLS
#Kegg
length(which(is.na(tempresultsKlevMA5$kegg_id)))
length(which(tempresultsMA5$KEGG=="None"))
length(which(!is.na(tempresultsKlevMA5$kegg_id)))
length(which(tempresultsMA5$KEGG!="None"))
#chebi
length(which(is.na(tempresultsKlevMA5$chebi_id)))
length(which(tempresultsMA5$ChEBI=="None"))
length(which(!is.na(tempresultsKlevMA5$chebi_id)))
length(which(tempresultsMA5$ChEBI!="None"))
#hmdb
length(which(is.na(tempresultsKlevMA5$hmdb_id)))
length(which(tempresultsMA5$HMDB=="None"))
length(which(!is.na(tempresultsKlevMA5$hmdb_id)))
length(which(tempresultsMA5$HMDB!="None"))
#pubchem
length(which(is.na(tempresultsKlevMA5$pubchem_id)))
length(which(tempresultsMA5$PubChem=="None"))
length(which(!is.na(tempresultsKlevMA5$pubchem_id)))
length(which(tempresultsMA5$PubChem!="None"))



#---------MetaboAnalyst5.0 and MetaFetcher using prostate cancer data--------------------------------------


query_results_table_cancer_MA5<-read.csv("MA5Comparison/MA_names_Priolo.csv")
rownames(query_results_table_cancer_MA5) <- query_results_table_cancer_MA5[,1]
query_results_table_cancer_MA5=as.data.frame(query_results_table_cancer_MA5)
NotFoundinMA5_cancer=which(!cancerData$Metabolite.%in% rownames(query_results_table_cancer_MA5))
NotFoundinMFMA5=which(!(rownames(query_results_table_cancer_MA5) %in% cancerData$Metabolite.[-NotFoundinMA5_cancer]))
tempresultsMF_cancer_MA5=cancerDataMFresults[-NotFoundinMA5_cancer,]
tempcancerDataMA5=cancerData[-NotFoundinMA5_cancer,]
tempresultsMA5_cancer=query_results_table_cancer_MA5[-NotFoundinMFMA5,]
#order based on names in klevIDs
tempresultsMA5_cancer=tempresultsMA5_cancer[match(as.character(tempcancerDataMA5$Metabolite.),as.character(rownames(tempresultsMA5_cancer))),]

tempresultsMA5_cancer=tempresultsMA5_cancer%>%
  mutate(across(everything(), as.character))
tempresultsMF_cancer_MA5=tempresultsMF_cancer_MA5 %>%
  mutate(across(everything(), as.character))
tempresultsMA5_cancer[is.na(tempresultsMA5_cancer)] <- "None"

mapply(compareIDs,tempresultsMA$kegg_id,tempKlevID$KEGG)

length(which(mapply(compareIDs,tempresultsMA5_cancer$KEGG,tempresultsMF_cancer_MA5$kegg_id)))

length(which(mapply(compareIDs,tempresultsMA5_cancer$HMDB,tempresultsMF_cancer_MA5$hmdb_id)))

length(which(mapply(compareIDs,tempresultsMA5_cancer$ChEBI,tempresultsMF_cancer_MA5$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))

length(which(mapply(compareIDs,tempresultsMA5_cancer$PubChem,tempresultsMF_cancer_MA5$pubchem_id)))


#NULLS
#Kegg
length(which(is.na(tempresultsMF_cancer_MA5$kegg_id)))
length(which(tempresultsMA5_cancer$KEGG=="None"))
length(which(!is.na(tempresultsMF_cancer_MA5$kegg_id)))
length(which(tempresultsMA5_cancer$KEGG!="None"))
#chebi
length(which(is.na(tempresultsMF_cancer_MA5$chebi_id)))
length(which(tempresultsMA5_cancer$ChEBI=="None"))
length(which(!is.na(tempresultsMF_cancer_MA5$chebi_id)))
length(which(tempresultsMA5_cancer$ChEBI!="None"))
#hmdb
length(which(is.na(tempresultsMF_cancer_MA5$hmdb_id)))
length(which(tempresultsMA5_cancer$HMDB=="None"))
length(which(!is.na(tempresultsMF_cancer_MA5$hmdb_id)))
length(which(tempresultsMA5_cancer$HMDB!="None"))
#pubchem
length(which(is.na(tempresultsMF_cancer_MA5$pubchem_id)))
length(which(tempresultsMA5_cancer$PubChem=="None"))
length(which(!is.na(tempresultsMF_cancer_MA5$pubchem_id)))
length(which(tempresultsMA5_cancer$PubChem!="None"))


