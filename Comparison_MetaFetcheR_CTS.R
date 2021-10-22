#Comparison CTS

cts_hmdb_input_MF<- read.xlsx("CTS_comparisons/CTSInputMetaFetcheR.xlsx",sheetName = "discoveryHMDB")
cts_hmdb_input_MF[c("kegg_id", "chebi_id","lipidmaps_id")] <- NA
cts_kegg_input_MF<- read.xlsx("CTS_comparisons/CTSInputMetaFetcheR.xlsx",sheetName = "discoveryKEGG")
cts_kegg_input_MF[c("hmdb_id", "chebi_id","lipidmaps_id")] <- NA
cts_lipidmaps_input_MF<- read.xlsx("CTS_comparisons/CTSInputMetaFetcheR.xlsx",sheetName = "discoveryLIPIDMAPS")
cts_lipidmaps_input_MF[c("hmdb_id", "chebi_id","kegg_id")] <- NA


#-------------------------------------------------------------------------------------------
cts_hmdb_output<- read.xlsx("CTS_comparisons/CTS-results.xlsx",sheetName = "D-HMDB")
cts_kegg_output<- read.xlsx("CTS_comparisons/CTS-results.xlsx",sheetName = "D-KEGG")
cts_lipidmaps_output<- read.xlsx("CTS_comparisons/CTS-results.xlsx",sheetName = "D-Lipidmaps")
cts_kegg_output_cancer<- read.xlsx("CTS_comparisons/CTS-results.xlsx",sheetName = "P-KEGG")
#For priolo et al we can use the already processed data for comparison
#cancerDataMFresults

cts_hmdb_output_MF<-resolve_metabolites(cts_hmdb_input_MF)

#unlist everything in the output of MF

cts_hmdb_output_MF<-transformMFoutput(cts_hmdb_output_MF$df,pubchem=FALSE)


cts_kegg_output_MF<-resolve_metabolites(cts_kegg_input_MF)

cts_kegg_output_MF<-transformMFoutput(cts_kegg_output_MF$df,pubchem=FALSE)


cts_lipidmaps_output_MF<-resolve_metabolites(cts_lipidmaps_input_MF)

cts_lipidmaps_output_MF<-transformMFoutput(cts_lipidmaps_output_MF$df,pubchem=FALSE)



#---------CTS and MetaFetcher using klev data (HMDB)---------------------------------------

#Matching and Non matching

mapply(compareIDs,tempresultsMA$kegg_id,tempKlevID$KEGG)
#removing the Chebi: from CTS and replacing NA with no result
tempChebi=as.factor(mapply(str_remove,cts_hmdb_output$ChEBI, "CHEBI:"))
tempChebi[which(is.na(tempChebi))]="No result"

length(which(mapply(compareIDs,cts_hmdb_output$KEGG,cts_hmdb_output_MF$kegg_id)))
length(which(mapply(compareIDs,cts_hmdb_output$LipidMAPS,cts_hmdb_output_MF$lipidmaps_id)))
length(which(mapply(compareIDs,tempChebi,cts_hmdb_output_MF$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))

#NULLS
#Kegg
length(which(is.na(cts_hmdb_output_MF$kegg_id)))
length(which(grepl("undified|No result",cts_hmdb_output$KEGG)))
length(which(!(is.na(cts_hmdb_output_MF$kegg_id))))
length(which(!(grepl("undified|No result",cts_hmdb_output$KEGG))))
#chebi
length(which(is.na(cts_hmdb_output_MF$chebi_id)))
length(which(grepl("undified|No result",tempChebi)))
length(which(!(is.na(cts_hmdb_output_MF$chebi_id))))
length(which(!(grepl("undified|No result",tempChebi))))

#lipidmaps
length(which(is.na(cts_hmdb_output_MF$lipidmaps_id)))
length(which(grepl("undified|No result",cts_hmdb_output$LipidMAPS)))
length(which(!is.na(cts_hmdb_output_MF$lipidmaps_id)))
length(which(!(grepl("undified|No result",cts_hmdb_output$LipidMAPS))))

#-------------------CTS and MetaFetcher using klev data (KEGG)---------------------------------------

#Matching and Non matching

#removing the Chebi: from CTS and replacing NA with no result
tempChebi=as.factor(mapply(str_remove,cts_kegg_output$ChEBI, "CHEBI:"))
tempChebi[which(is.na(tempChebi))]="No result"

length(which(mapply(compareIDs,cts_kegg_output$Human.Metabolome.Database,cts_kegg_output_MF$hmdb_id)))
length(which(mapply(compareIDs,cts_kegg_output$LipidMAPS,cts_kegg_output_MF$lipidmaps_id)))
length(which(mapply(compareIDs,tempChebi,cts_kegg_output_MF$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))

#NULLS
#HMDB
length(which(is.na(cts_kegg_output_MF$hmdb_id)))
length(which(grepl("undified|No result",cts_kegg_output$Human.Metabolome.Database)))
length(which(!(is.na(cts_kegg_output_MF$hmdb_id))))
length(which(!(grepl("undified|No result",cts_kegg_output$Human.Metabolome.Database))))
#lipidmaps
length(which(is.na(cts_kegg_output_MF$lipidmaps_id)))
length(which(grepl("undified|No result",cts_kegg_output$LipidMAPS)))
length(which(!(is.na(cts_kegg_output_MF$lipidmaps_id))))
length(which(!(grepl("undified|No result",cts_kegg_output$LipidMAPS))))

#chebi
length(which(is.na(cts_kegg_output_MF$chebi_id)))
length(which(grepl("undified|No result",tempChebi)))
length(which(!(is.na(cts_kegg_output_MF$chebi_id))))
length(which(!(grepl("undified|No result",tempChebi))))





#-------------------CTS and MetaFetcher using klev data (LIPIDMAPS)---------------------------------------

#Matching and Non matching

#removing the Chebi: from CTS and replacing NA with no result
tempChebi=as.factor(mapply(str_remove,cts_lipidmaps_output$ChEBI, "CHEBI:"))
tempChebi[which(is.na(tempChebi))]="No result"

length(which(mapply(compareIDs,cts_lipidmaps_output$Human.Metabolome.Database,cts_lipidmaps_output_MF$hmdb_id)))
length(which(mapply(compareIDs,cts_lipidmaps_output$KEGG,cts_lipidmaps_output_MF$kegg_id)))
length(which(mapply(compareIDs,tempChebi,cts_lipidmaps_output_MF$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))

#NULLS
#HMDB
length(which(is.na(cts_lipidmaps_output_MF$hmdb_id)))
length(which(grepl("undified|No result",cts_lipidmaps_output$Human.Metabolome.Database)))
length(which(!(is.na(cts_lipidmaps_output_MF$hmdb_id))))
length(which(!(grepl("undified|No result",cts_lipidmaps_output$Human.Metabolome.Database))))
#kegg
length(which(is.na(cts_lipidmaps_output_MF$kegg_id)))
length(which(grepl("undified|No result",cts_lipidmaps_output$KEGG)))
length(which(!(is.na(cts_lipidmaps_output_MF$kegg_id))))
length(which(!(grepl("undified|No result",cts_lipidmaps_output$KEGG))))

#chebi
length(which(is.na(cts_lipidmaps_output_MF$chebi_id)))
length(which(grepl("undified|No result",tempChebi)))
length(which(!(is.na(cts_lipidmaps_output_MF$chebi_id))))
length(which(!(grepl("undified|No result",tempChebi))))



#-----------------Comparison 2 CTS and Metafetcher (Cancer data)------------------------------
#cancerDataMFresults
#cts_kegg_output_cancer
#Matching and Non matching


mapply(compareIDs,tempresultsMA$kegg_id,tempKlevID$KEGG)
#removing the Chebi: from CTS and replacing NA with no result
tempChebi=as.factor(mapply(str_remove,cts_kegg_output_cancer$ChEBI, "CHEBI:"))
tempChebi[which(is.na(tempChebi))]="No result"
#Matching and Non matching


length(which(mapply(compareIDs,cts_kegg_output_cancer$Human.Metabolome.Database,cancerDataMFresults$hmdb_id)))
length(which(mapply(compareIDs,cts_kegg_output_cancer$LipidMAPS,cancerDataMFresults$lipidmaps_id)))
length(which(mapply(compareIDs,tempChebi,cancerDataMFresults$chebi_id)))
#length(which(mapply(compareIDs,sstempresultsMA$lipidmaps_id,tempresultsKlev$lipidmaps_id)))

#NULLS
#HMDB
length(which(is.na(cancerDataMFresults$hmdb_id)))
length(which(grepl("undified|No result",cts_kegg_output_cancer$Human.Metabolome.Database)))
length(which(!(is.na(cancerDataMFresults$hmdb_id))))
length(which(!(grepl("undified|No result",cts_kegg_output_cancer$Human.Metabolome.Database))))
#lipidmaps
length(which(is.na(cancerDataMFresults$lipidmaps_id)))
length(which(grepl("undified|No result",cts_kegg_output_cancer$LipidMAPS)))
length(which(!(is.na(cancerDataMFresults$lipidmaps_id))))
length(which(!(grepl("undified|No result",cts_kegg_output_cancer$LipidMAPS))))

#chebi
length(which(is.na(cancerDataMFresults$chebi_id)))
length(which(grepl("undified|No result",tempChebi)))
length(which(!(is.na(cancerDataMFresults$chebi_id))))
length(which(!(grepl("undified|No result",tempChebi))))




