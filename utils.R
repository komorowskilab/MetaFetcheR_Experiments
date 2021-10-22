compareIDs<-function(var1,var2)
{
  
  if (var1!="-")
    grepl(paste("?",var1,"?",sep=""),var2)
  else
    return(FALSE)
}
transformMFoutput<-function(output,hmdb=TRUE,chebi=TRUE,pubchem=TRUE,lipidmaps=TRUE,kegg=TRUE)
{
  temp=NULL
  if(kegg==TRUE)
    temp$kegg_id=list_to_text(output$kegg_id,sep="-")
  if(hmdb==TRUE)
    temp$hmdb_id=list_to_text(output$hmdb_id,sep="-")
  if(chebi==TRUE)
    temp$chebi_id=list_to_text(output$chebi_id,sep="-")
  if(pubchem==TRUE)
    temp$pubchem_id=list_to_text(output$pubchem_id,sep="-")
  if(lipidmaps==TRUE)
    temp$lipidmaps_id=list_to_text(output$lipidmaps_id,sep="-")
  temp=as.data.frame(temp)
  return(temp)
  
}