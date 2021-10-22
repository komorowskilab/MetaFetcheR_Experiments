library(MetaboAnalystR)
library(metafetcher)
library(xlsx)
library(dplyr)
library(ggplot2)
#---------------for Comparing MetaFetcheR with MetaboAnalyst results using Diamanti et al (case 2) and also (Metafetcher results and manual curation)(case 1) -------------------------
#This following script is used for for Comparison 1 Comparing with MetaboAnalyst results using Diamanti et al (case 2) and also (Metafetcher results and manual curation)(case 1)

KlevIds=read.xlsx("20190426_metabolite_ids.xlsx",sheetIndex = 1)
KlevIDSSelected=KlevIds[,c("KEGG","HMDB","ChEBI","PubChem","LIPIDMAPS")]
df_from_csv <- read.csv("discovery.csv", stringsAsFactors=FALSE)
colnames(KlevIDSSelected)=c("kegg_id","hmdb_id","chebi_id","pubchem_id","lipidmaps_id")

resultsKlev=read.csv("results.csv",sep = ";")

name.vec<-paste0(trimws(KlevIds$common.name), collapse=";")
toSend = list(queryList = name.vec, inputType = "name")
call="http://api.xialab.ca/mapcompounds"
query_results<-httr::POST(call, body = toSend, encode = "json")
query_results_text <- content(query_results, "text", encoding = "UTF-8")
query_results_json <-RJSONIO::fromJSON(query_results_text, flatten = TRUE)
query_results_table <- t(rbind.data.frame(query_results_json))
query_results_tableMA5<-read.csv("MA5Comparison/MA_names.csv")
rownames(query_results_table) <- query_results_table[,1]
query_results_table=as.data.frame(query_results_table)
NotFoundinMA=which(!(trimws(KlevIds$common.name) %in% rownames(query_results_table)))
NotFoundinKlev=which(!(rownames(query_results_table) %in% trimws(KlevIds$common.name)[-NotFoundinMA]))
tempresultsKlev=resultsKlev[-NotFoundinMA,]
tempKlevID=KlevIds[-NotFoundinMA,]
tempresultsMA=query_results_table[-NotFoundinKlev,]
#order based on names in klevIDs
tempresultsMA=tempresultsMA[match(tempKlevID$common.name,rownames(tempresultsMA)),]
#---------------------------------------------------------------------
#for Comparing  MetaboAnalystresults and MetafetcheR using Prostate cancer data (case 2)
#---------------------------------------------------------------------
cancerData <- read.xlsx("prostate_cancer_dataset.xlsx",sheetName = "Human Tumors")
cancerDataMF=NULL
cancerDataMF$kegg_id=cancerData$KEGG.ID
cancerDataMF$hmdb_id=rep(NA,length(cancerData$KEGG.ID))
cancerDataMF$chebi_id=rep(NA,length(cancerData$KEGG.ID))
cancerDataMF$pubchem_id=rep(NA,length(cancerData$KEGG.ID))
cancerDataMF$lipidmaps_id=rep(NA,length(cancerData$KEGG.ID))
cancerDataMF=apply(cancerDataMF,2,as.character)
cancerDataMF=as.data.frame(cancerDataMF)
cancerDataMFresults <- resolve_metabolites(cancerDataMF)
temp=NULL
temp$kegg_id=list_to_text(cancerDataMFresults$df$kegg_id,sep="-")
temp$hmdb_id=list_to_text(cancerDataMFresults$df$hmdb_id,sep="-")
temp$chebi_id=list_to_text(cancerDataMFresults$df$chebi_id,sep="-")
temp$pubchem_id=list_to_text(cancerDataMFresults$df$pubchem_id,sep="-")
temp$lipidmaps_id=list_to_text(cancerDataMFresults$df$lipidmaps_id,sep="-")
temp=as.data.frame(temp)
cancerDataMFresults=temp

colnames(cancerDataMFresults)=c("kegg_id","hmdb_id","chebi_id","pubchem_id","lipidmaps_id")

name.vec<-paste0(trimws(cancerData$Metabolite.), collapse=";")
toSend = list(queryList = name.vec, inputType = "name")
call="http://api.xialab.ca/mapcompounds"
query_results<-httr::POST(call, body = toSend, encode = "json")
query_results_text <- content(query_results, "text", encoding = "UTF-8")
query_results_json <- RJSONIO::fromJSON(query_results_text, flatten = TRUE)
query_results_table_cancer <- t(rbind.data.frame(query_results_json))
rownames(query_results_table_cancer) <- query_results_table_cancer[,1]
query_results_table_cancer=as.data.frame(query_results_table_cancer)
NotFoundinMA_cancer=which(!cancerData$Metabolite.%in% rownames(query_results_table_cancer))
NotFoundinMF=which(!(rownames(query_results_table_cancer) %in% cancerData$Metabolite.[-NotFoundinMA_cancer]))
tempresultsMF_cancer=cancerDataMFresults[-NotFoundinMA_cancer,]
tempcancerData=cancerData[-NotFoundinMA_cancer,]
tempresultsMA_cancer=query_results_table_cancer[-NotFoundinMF,]
#order based on names in klevIDs
tempresultsMA_cancer=tempresultsMA_cancer[match(as.character(tempcancerData$Metabolite.),as.character(rownames(tempresultsMA_cancer))),]


mapply(compareIDs,tempresultsMA_cancer$kegg_id,tempresultsMF_cancer$kegg_id)
mapply(compareIDs,tempresultsMA$kegg_id,tempKlevID$KEGG)

#---------MS_Targeted+manual curation and MetaFetcher using Diamanti et al data(case 1)---------------------------------------
#NULLS
#Kegg
length(which(is.na(resultsKlev$kegg_id)))
length(which(resultsKlev$fKEGG==""))
length(which(!resultsKlev$fKEGG==""))
length(which(is.na(KlevIds$KEGG)))
#chebi
length(which(is.na(resultsKlev$chebi_id)))
length(which(is.na(resultsKlev$fChEBI)))
length(which(!is.na(resultsKlev$fChEBI)))
length(which(is.na(KlevIds$ChEBI)))
#hmdb
length(which(is.na(resultsKlev$hmdb_id)))
length(which(is.na(resultsKlev$fHMDB)))
length(which(!is.na(resultsKlev$fHMDB)))
length(which(is.na(KlevIds$HMDB)))
#lipidmaps
length(which(is.na(resultsKlev$lipidmaps_id)))
length(which(is.na(resultsKlev$fLIPIDMAPS)))
length(which(!is.na(resultsKlev$fLIPIDMAPS)))
length(which(is.na(KlevIds$LIPIDMAPS)))
#pubchem

length(which(is.na(resultsKlev$fPubChem)))
length(which(!is.na(resultsKlev$fPubChem)))
length(which(is.na(KlevIds$PubChem)))
#comparison MS_targeted
length(which(mapply(compareIDs,resultsKlev$fKEGG,resultsKlev$kegg_id)))
length(which(mapply(compareIDs,resultsKlev$fHMDB,resultsKlev$hmdb_id)))
length(which(mapply(compareIDs,resultsKlev$fChEBI,resultsKlev$chebi_id)))
length(which(mapply(compareIDs,resultsKlev$fLIPIDMAPS,resultsKlev$lipidmaps_id)))
length(which(mapply(compareIDs,resultsKlev$fPubChem,resultsKlev$pubchem_id)))



S1=read.xlsx("Source1.xlsx",sheetIndex = 1)
rownames(S1)=S1[,1]
S1=S1[,-1]

#---------MetaboAnalyst and MetaFetcher using Diamanti et al data (case 2)---------------------------------------


mapply(compareIDs,tempresultsMA$kegg_id,tempKlevID$KEGG)

length(which(mapply(compareIDs,tempresultsMA$kegg_id,tempresultsKlev$kegg_id)))
length(which(mapply(compareIDs,tempresultsMA$hmdb_id,tempresultsKlev$hmdb_id)))
length(which(mapply(compareIDs,tempresultsMA$chebi_id,tempresultsKlev$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))
length(which(mapply(compareIDs,tempresultsMA$pubchem_id,tempresultsKlev$pubchem_id)))

#----MetaboAnalyst------------------
#NULLS
#Kegg
length(which(is.na(tempresultsKlev$kegg_id)))
length(which(tempresultsMA$kegg_id=="-"))
length(which(!is.na(tempresultsKlev$kegg_id)))
length(which(tempresultsMA$kegg_id!="-"))
#chebi
length(which(is.na(tempresultsKlev$chebi_id)))
length(which(tempresultsMA$chebi_id=="-"))
length(which(!is.na(tempresultsKlev$chebi_id)))
length(which(tempresultsMA$chebi_id!="-"))
#hmdb
length(which(is.na(tempresultsKlev$hmdb_id)))
length(which(tempresultsMA$hmdb_id=="-"))
length(which(!is.na(tempresultsKlev$hmdb_id)))
length(which(tempresultsMA$hmdb_id!="-"))
#pubchem
length(which(is.na(tempresultsKlev$pubchem_id)))
length(which(tempresultsMA$pubchem_id=="-"))
length(which(!is.na(tempresultsKlev$pubchem_id)))
length(which(tempresultsMA$pubchem_id!="-"))


S2=read.xlsx("Source1.xlsx",sheetIndex = 4)
rownames(S2)=S2[,1]
S2=S2[,-1]

#---------MetaboAnalyst and MetaFetcher using prostate cancer data (case 2)--------------------------------------
mapply(compareIDs,tempresultsMA$kegg_id,tempKlevID$KEGG)

length(which(mapply(compareIDs,tempresultsMA_cancer$kegg_id,tempresultsMF_cancer$kegg_id)))

length(which(mapply(compareIDs,tempresultsMA_cancer$hmdb_id,tempresultsMF_cancer$hmdb_id)))

length(which(mapply(compareIDs,tempresultsMA_cancer$chebi_id,tempresultsMF_cancer$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))

length(which(mapply(compareIDs,tempresultsMA_cancer$pubchem_id,tempresultsMF_cancer$pubchem_id)))

#NULLS
#Kegg
length(which(is.na(tempresultsMF_cancer$kegg_id)))
length(which(tempresultsMA_cancer$kegg_id=="-"))
length(which(!is.na(tempresultsMF_cancer$kegg_id)))
length(which(tempresultsMA_cancer$kegg_id!="-"))
#chebi
length(which(is.na(tempresultsMF_cancer$chebi_id)))
length(which(tempresultsMA_cancer$chebi_id=="-"))
length(which(!is.na(tempresultsMF_cancer$chebi_id)))
length(which(tempresultsMA_cancer$chebi_id!="-"))
#hmdb
length(which(is.na(tempresultsMF_cancer$hmdb_id)))
length(which(tempresultsMA_cancer$hmdb_id=="-"))
length(which(!is.na(tempresultsMF_cancer$hmdb_id)))
length(which(tempresultsMA_cancer$hmdb_id!="-"))
#pubchem
length(which(is.na(tempresultsMF_cancer$pubchem_id)))
length(which(tempresultsMA_cancer$pubchem_id=="-"))
length(which(!is.na(tempresultsMF_cancer$pubchem_id)))
length(which(tempresultsMA_cancer$pubchem_id!="-"))


S3=read.xlsx("Source1.xlsx",sheetIndex = 5)
rownames(S3)=S3[,1]
S3=S3[,-1]

