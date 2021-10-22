library(ggsci)
library(wesanderson)

#----------------------------Figure 4-------------------------------------------------------------------------------
#We used function do_consistency_test in the R package but we didn't export it . We will make it public in the MetaFetcheR package.You can load the enviroment ResultsCoverage.RData that we provide in the repository
resulthmdb2=as.data.frame(resulthmdb2)
resultchebi2=as.data.frame(resultchebi2)
resultpubchem2=as.data.frame(resultpubchem2)
resultkegg2=as.data.frame(resultkegg2)
resultlipidmaps2=as.data.frame(resultlipidmaps2)
x=rep(1:100)
decision=append(append(rep("HMDB",100),rep("ChEBI",100)),append(rep("LIPID MAPS",100),append(rep("PubChem",100),rep("KEGG",100))))

AllResolvedAtt=as.data.frame(cbind(resulthmdb2$`Resolved attributes`,resultchebi2$`Resolved attributes`,resultlipidmaps2$`Resolved attributes`,resultpubchem2$`Resolved attributes`,resultkegg2$`Resolved attributes`))
colnames(AllResolvedAtt)=c("HMDB","ChEBI","LIPID MAPS","PubChem","KEGG")

AllAmbAtt=as.data.frame(cbind(resulthmdb2$`Ambigous attributes`,resultchebi2$`Ambigous attributes`,resultlipidmaps2$`Ambigous attributes`,resultpubchem2$`Ambigous attributes`,resultkegg2$`Ambigous attributes`))
colnames(AllAmbAtt)=c("HMDB","ChEBI","LIPID MAPS","PubChem","KEGG")


AllUnresolvedAtt=as.data.frame(cbind(resulthmdb2$`Missing attributes`,resultchebi2$`Missing attributes`,resultlipidmaps2$`Missing attributes`,resultpubchem2$`Missing attributes`,resultkegg2$`Missing attributes`))
colnames(AllUnresolvedAtt)=c("HMDB","ChEBI","LIPID MAPS","PubChem","KEGG")

ggplot(AllResolvedAtt,aes(x=iterations)) +
  geom_line(aes(y = hmdb), color = "darkorange") +
  geom_line(aes(y =chebi), color="darkblue")+
  geom_line(aes(y =lipidmaps), color="purple")+
  geom_line(aes(y =pubchem), color="darkgreen")+
  geom_line(aes(y =kegg), color="darkred")


df1.m=reshape2::melt(AllResolvedAtt)
df1.m$value=as.numeric(df1.m$value)
p1<-ggplot(as.data.frame(df1.m) ,aes(x=variable, y=value, fill=decision)) +
        geom_boxplot()+ylab("percentage")+xlab("")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                             panel.grid.minor = element_blank(),axis.title = element_text(size = 20),axis.text.y=element_text(size=20),axis.text.x=element_text(angle = 45, hjust = 1,size=20),legend.text=element_text(size=20),axis.line = element_line(colour = "grey"),panel.background = element_blank())+scale_fill_npg()

df2.m=reshape2::melt(AllAmbAtt)
df2.m$value=as.numeric(df2.m$value)
p2<-ggplot(as.data.frame(df2.m) ,aes(x=variable, y=value, fill=decision)) +
        geom_boxplot()+ylab("percentage")+xlab("")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                             panel.grid.minor = element_blank(),axis.title = element_text(size = 20),axis.text.y=element_text(size=20),axis.text.x=element_text(angle = 45, hjust = 1,size=20),legend.text=element_text(size=20),axis.line = element_line(colour = "grey"),panel.background = element_blank())+scale_fill_npg()



df3.m=reshape2::melt(AllUnresolvedAtt)
df3.m$value=as.numeric(df3.m$value)
p3<-ggplot(as.data.frame(df3.m) ,aes(x=variable, y=value, fill=decision)) +
        geom_boxplot()+ylab("percentage")+xlab("")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                             panel.grid.minor = element_blank(),axis.title = element_text(size = 20),axis.text.y=element_text(size=20),axis.text.x=element_text(angle = 45, hjust = 1,size=20),legend.text=element_text(size=20),axis.line = element_line(colour = "grey"),panel.background = element_blank())+scale_fill_npg()


#cowplot::plot_grid(p1, p2, p3,labels = "AUTO")
cowplot::plot_grid(p1, p2,p3,labels = "AUTO", ncol = 3)
prow <- cowplot::plot_grid(
  p1 + theme(legend.position="none"),
  p2 + theme(legend.position="none"),
  p3 + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "B", "C"),
  hjust = -1,
  nrow = 1
)
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of
# the width of one plot (via rel_widths).
cowplot::plot_grid(prow, legend, rel_widths = c(4, .4))

#--------------------Comparison 1-MS_Targeted with klev data (Supplementary Figure 6)-------------------------
pal <- wes_palette("GrandBudapest1", 100, type = "continuous")

S1=read.xlsx("Source1.xlsx",sheetIndex = 1)
rownames(S1)=S1[,1]
S1=S1[,-1]
S1=S1[,c(1,3,5,6)]
colnames(S1)=c("Matching","MS_Targeted/Manually curated","Empty feilds in MS_Targeted","Not Matching")
S1=S1[,c(1,2)]
S1melt=reshape2::melt(S1)
S1melt$DB=rep( c("HMDB","KEGG","ChEBI","LIPID MAPS","PubChem"),2)
S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated"))
#S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
#S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
S1melt$DB<-factor(S1melt$DB, levels = c("HMDB","KEGG","ChEBI","LIPID MAPS","PubChem"))
S1melt$variable=as.character(S1melt$variable)

ggplot(data=S1melt, aes(x=DB,y=value,fill=variable)) +
  geom_bar(stat="identity")+ylab("")+xlab("")+guides(fill=guide_legend(title=""))+theme(axis.ticks = element_blank(),text = element_text(size=14), axis.text.y = element_blank())+scale_fill_npg()
  #scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))


#--------------------Comparison 2-MetaboAnalyst with klev data-------------------------

S2=read.xlsx("Source1.xlsx",sheetIndex = 4)
rownames(S2)=S2[,1]
S2=S2[,-1]
colnames(S2)=c("Matching","Empty MetafetcheR","Empty MetaboAnalyst","Non empty MetafetcheR","Non empty MetaboAnalyst","Not Matching")

S2melt=reshape2::melt(S2,id.vars="id")
S2melt$DB=rep( c("HMDB","KEGG","ChEBI","PubChem"),6)
S2melt$variable <- factor(S2melt$variable, levels = c("Empty MetafetcheR","Empty MetaboAnalyst","Non empty MetafetcheR","Non empty MetaboAnalyst","Matching","Not Matching"))
#S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
#S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
S2melt$DB<-factor(S2melt$DB, levels =  c("HMDB","KEGG","ChEBI","LIPID MAPS","PubChem"))
S2melt$variable=as.character(S2melt$variable)

ggplot(data=S2melt, aes(x=DB,y=value,fill=variable)) +
  geom_bar(stat="identity")+ylab("")+xlab("")+guides(fill=guide_legend(title=""))+theme(axis.ticks = element_blank(), axis.text.y = element_blank())+scale_fill_npg()

  #scale_fill_brewer(palette = "Blues")
ggdotchart(S2melt, x = "DB", y = "value",
           color = "variable",                                # Color by groups
          # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "variable",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(S2melt$value,1),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
)
  geom_hline(yintercept = 0, linetype = 2, color = "lightgray")


#--------------------Comparison 1-MS_Targeted with klev data (Supplementary Figure 6)-------------------------
S1=read.xlsx("Source1.xlsx",sheetIndex = 1)
rownames(S1)=S1[,1]
S1=S1[,-1]
S1=S1[,c(1,3,5,6)]
colnames(S1)=c("Matching","MS_targeted/Manually curated","Empty feilds in MS_Targeted","Not Matching")
S1=S1[,c(1,2)]
S1$id=c("HMDB","KEGG","ChEBI","LIPID MAPS","PubChem")
S1melt=reshape2::melt(S1,id.vars="id")
colnames(S1melt)<-c("DB","variable","value")
#S1melt$DB=rep( c("HMDB","KEGG","ChEBI","LIPID MAPS","PubChem"),2)
S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_targeted/Manually curated"))
#S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
#S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
S1melt$DB<-factor(S1melt$DB, levels = c("PubChem","LIPID MAPS","ChEBI","KEGG","HMDB"))
S1melt$variable=as.character(S1melt$variable)

ggplot(data=S1melt, aes(x=DB,y=value,fill=variable)) +
  geom_bar(stat="identity")+ylab("")+xlab("")+guides(fill=guide_legend(title=""))+theme(axis.ticks = element_blank(),text = element_text(size=14), axis.text.y = element_blank())+scale_fill_npg()
#scale_fill_manual(values = wes_palette("GrandBudapest1", n = 3))

pS3 <- ggplot(S1melt, aes(DB, value,fill=variable))
pgS3 <- pS3+ geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() + xlab("")+ylab("")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),panel.background = element_blank())+scale_fill_npg()

print(pgS3)


#----------Case 2-MetaboAnalystR with Diamanti et al data (Figure.2 A)---------------
S2=read.xlsx("Source1.xlsx",sheetIndex = 4)
rownames(S2)=S2[,1]
S2=S2[,-1]
colnames(S2)=c("Matching","Empty in MetaFetcheR","Empty in MetaboAnalyst","Non empty in MetaFetcheR","Non empty in MetaboAnalyst","Not Matching")
S2$id=c("HMDB","KEGG","ChEBI","PubChem")
S2melt=reshape2::melt(S2,id.vars="id")
colnames(S2melt)<-c("DB","variable","value")


#S2melt$Var2=rep( c("HMDB","KEGG","ChEBI","PubChem"),6)
S2melt$variable <- factor(S2melt$variable, levels = c("Empty in MetaFetcheR","Empty in MetaboAnalyst","Non empty in MetaFetcheR","Non empty in MetaboAnalyst","Matching","Not Matching"))
#S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
#S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
#S2melt$DB<-factor(S2melt$DB, levels =  c("HMDB","KEGG","ChEBI","LIPID MAPS","PubChem"))
S2melt$DB<-factor(S2melt$DB, levels =c("PubChem","LIPID MAPS","ChEBI","KEGG","HMDB"))
#S2melt$variable=as.character(S2melt$variable)


ggdotchart(S2melt, x = "DB", y = "value",
           color = "variable",                                # Color by groups
           # palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           group = "variable",                                # Order by groups
           dot.size = 6,                                 # Large dot size
           label = round(S2melt$value,1),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
)
geom_hline(yintercept = 0, linetype = 2, color = "lightgray")




#Lattice package trial
 barchart(DB ~ value | factor(variable), data=S2melt,
           main="barchart",
           scales=list(cex=0.5),
           layout=c(3, 1))

 #Doing the same thing as lattice but with ggplot

 p <- ggplot(S2melt, aes(DB, value,fill=variable))
 pg <- p + geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() +ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
 panel.grid.minor = element_blank(),text = element_text(size=14),axis.line = element_line(colour = "grey"),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

print(pg)

#----------case 2-MetaboAnalystR with Priolo et al data (Figure 2 C)--------------------------------

    S3=read.xlsx("Source1.xlsx",sheetIndex = 5)
    rownames(S3)=S3[,1]
    S3=S3[,-1]
    colnames(S3)=c("Matching","Not Matching","Empty in MetaFetcheR","Empty in MetaboAnalyst","Non empty in MetaFetcheR","Non empty in MetaboAnalyst")
    S3$id=c("HMDB","KEGG","ChEBI","PubChem")
    S3melt=reshape2::melt(S3,id.vars="id")
    colnames(S3melt)<-c("DB","variable","value")


    S3melt$variable <- factor(S3melt$variable, levels = c("Empty in MetaFetcheR","Empty in MetaboAnalyst","Non empty in MetaFetcheR","Non empty in MetaboAnalyst","Matching","Not Matching"))
    #S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
    #S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
    S3melt$DB<-factor(S3melt$DB, levels =  c("PubChem","LIPID MAPS","ChEBI","KEGG","HMDB"))
    S3melt$variable=as.character(S3melt$variable)

    #Doing the same thing as lattice but with ggplot

    pS5 <- ggplot(S3melt, aes(DB, value,fill=variable))
    pgS5 <- pS5 + geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() +ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                        panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

    print(pgS5)


 #----------Comparison 4a-CTS with MetaFetcheR using Diamanti et al data--Starting with HMDB (Supplementary Figure S4-A)---------------

    S4a=read.xlsx("Source2.xlsx",sheetIndex = 16)
    S4a=S4a[1:6,1:4]
    rownames(S4a)=S4a[,1]
    S4a=S4a[,-1]
    S4a=t(S4a)
    colnames(S4a)=c("Matching","Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Not Matching")
    S4a=as.data.frame(S4a)
    S4a$id=c("KEGG","ChEBI","LIPID MAPS")
    S4amelt=reshape2::melt(S4a,id.vars="id")
    colnames(S4amelt)<-c("DB","variable","value")


    S4amelt$variable <- factor(S4amelt$variable, levels = c("Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Matching","Not Matching"))
    #S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
    #S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
    S4amelt$DB<-factor(S4amelt$DB, levels =  c("LIPID MAPS","ChEBI","KEGG"))
   # S4amelt$variable=as.character(S4amelt$variable)

    #Doing the same thing as lattice but with ggplot

    pS6a <- ggplot(S4amelt, aes(DB, value,fill=variable))
    pgS6a <- pS6a + geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() + ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                            panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

    print(pgS6a)

    #----------Comparison 4b-CTS with MetaFetcheR using Diamanti et al data different plot)----Starting with KEGG(Supplementary Figure S4-B)-------------------------------------------

    S4b=read.xlsx("Source2.xlsx",sheetIndex = 17)
    S4b=S4b[1:6,1:4]
    rownames(S4b)=S4b[,1]
    S4b=S4b[,-1]
    S4b=t(S4b)
    colnames(S4b)=c("Matching","Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Not Matching")
    S4b=as.data.frame(S4b)
    S4b$id=c("HMDB","ChEBI","LIPID MAPS")
    S4bmelt=reshape2::melt(S4b,id.vars="id")
    colnames(S4bmelt)<-c("DB","variable","value")


    S4bmelt$variable <- factor(S4bmelt$variable, levels = c("Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Matching","Not Matching"))
    #S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
    #S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
    S4bmelt$DB<-factor(S4bmelt$DB, levels =  c("LIPID MAPS","ChEBI","HMDB"))
   # S4amelt$variable=as.character(S4melt$variable)

    #Doing the same thing as lattice but with ggplot

    pS6b <- ggplot(S4bmelt, aes(DB, value,fill=variable))
    pgS6b<- pS6b + geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() + ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                              panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

    print(pgS6b)
    #----------Comparison 4c-CTS with MetaFetcheR using Diamanti et al data different plot)----Starting with LIPID MAPS(Supplementary Figure S4-c)-----------

    S4c=read.xlsx("Source2.xlsx",sheetIndex = 18)
    S4c=S4c[1:6,1:4]
    rownames(S4c)=S4c[,1]
    S4c=S4c[,-1]
    S4c=t(S4c)
    colnames(S4c)=c("Matching","Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Not Matching")
    S4c=as.data.frame(S4c)
    S4c$id=c("HMDB","KEGG","ChEBI")
    S4cmelt=reshape2::melt(S4c,id.vars="id")
    colnames(S4cmelt)<-c("DB","variable","value")


    S4cmelt$variable <- factor(S4cmelt$variable, levels = c("Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Matching","Not Matching"))
    #S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
    #S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
    S4cmelt$DB<-factor(S4cmelt$DB, levels =  c("HMDB","ChEBI","KEGG"))
    # S4amelt$variable=as.character(S4melt$variable)

    #Doing the same thing as lattice but with ggplot

    pS6c <- ggplot(S4cmelt, aes(DB, value,fill=variable))
    pgS6c<- pS6c + geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() + ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                             panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

    print(pgS6c)

    #----------Comparison 5-CTS with MetaFetcheR using Priolo et al data different plot)----Starting with LIPID MAPS-(Supplementary Figure S5)-----------

    S5=read.xlsx("Source2.xlsx",sheetIndex = 19)
    S5=S5[1:6,1:4]
    rownames(S5)=S5[,1]
    S5=S5[,-1]
    S5=t(S5)
    colnames(S5)=c("Matching","Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Not Matching")
    S5=as.data.frame(S5)
    S5$id=c("HMDB","ChEBI","LIPID MAPS")
    S5melt=reshape2::melt(S5,id.vars="id")
    colnames(S5melt)<-c("DB","variable","value")


    S5melt$variable <- factor(S5melt$variable, levels = c("Empty in MetaFetcheR","Empty in CTS","Non empty in MetaFetcheR","Non empty in CTS","Matching","Not Matching"))
    #S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
    #S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
    S5melt$DB<-factor(S5melt$DB, levels =  c("LIPID MAPS","ChEBI","HMDB"))
    # S4amelt$variable=as.character(S4melt$variable)

    #Doing the same thing as lattice but with ggplot

    pS7 <- ggplot(S5melt, aes(DB, value,fill=variable))
    pgS7<-   pS7 + geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() + ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                                                                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

    print(pgS7)
    #-------------------------------Comparison 6a-MA5 with MetaFetcheR using Diamanti et al)(Figure 2 C)-------------------------

    S6=read.xlsx("Source2.xlsx",sheetIndex = 12)
    S6=S6[1:6,1:5]
    rownames(S6)=S6[,1]
    S6=S6[,-1]
    S6=t(S6)
    colnames(S6)=c("Matching","Empty in MetaFetcheR","Empty in MetaboAnalystWeb","Non empty in MetaFetcheR","Non empty in MetaboAnalystWeb","Not Matching")
    S6=as.data.frame(S6)
    S6$id=c("HMDB","KEGG","ChEBI","PubChem")
    S6melt=reshape2::melt(S6,id.vars="id")
    colnames(S6melt)<-c("DB","variable","value")


    S6melt$variable <- factor(S6melt$variable, levels = c("Empty in MetaFetcheR","Empty in MetaboAnalystWeb","Non empty in MetaFetcheR","Non empty in MetaboAnalystWeb","Matching","Not Matching"))
    #S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
    #S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
   # S6melt$DB<-factor(S6melt$DB, levels = c("HMDB","KEGG","ChEBI","PubChem"))
    S6melt$DB<-factor(S6melt$DB, levels = c("PubChem","KEGG","ChEBI","HMDB"))
    # S4amelt$variable=as.character(S4melt$variable)

    #Doing the same thing as lattice but with ggplot

    pS8 <- ggplot(S6melt, aes(DB, value,fill=variable))
    pgS8<-   pS8+ geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() + ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                                                                    panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

    print(pgS8)


   #---------------------------------------Comparison 6b-MA5 with MetaFetcheR using Priolo et al data (Figure 2 D)----------------------------

    S7=read.xlsx("Source2.xlsx",sheetIndex = 13)
    S7=S7[1:6,1:5]
    rownames(S7)=S7[,1]
    S7=S7[,-1]
    S7=t(S7)
    colnames(S7)=c("Matching","Empty in MetaFetcheR","Empty in MetaboAnalystWeb","Non empty in MetaFetcheR","Non empty in MetaboAnalystWeb","Not Matching")
    S7=as.data.frame(S7)
    S7$id=c("HMDB","KEGG","ChEBI","PubChem")
    S7melt=reshape2::melt(S7,id.vars="id")
    colnames(S7melt)<-c("DB","variable","value")


    S7melt$variable <- factor(S7melt$variable, levels = c("Empty in MetaFetcheR","Empty in MetaboAnalystWeb","Non empty in MetaFetcheR","Non empty in MetaboAnalystWeb","Matching","Not Matching"))
    #S1melt$DB=rep(c("hmdb","kegg","chebi","lipidmaps","pubchem"),4)
    #S1melt$variable <- factor(S1melt$variable, levels = c("Matching","MS_Targeted/Manually curated","Empty feids MS_Targeted","Not Matching"))
    # S6melt$DB<-factor(S6melt$DB, levels = c("HMDB","KEGG","ChEBI","PubChem"))
    S7melt$DB<-factor(S7melt$DB, levels = c("PubChem","KEGG","ChEBI","HMDB"))
    # S4amelt$variable=as.character(S4melt$variable)

    #Doing the same thing as lattice but with ggplot

    pS9 <- ggplot(S7melt, aes(DB, value,fill=variable))
    pgS9<-   pS9+ geom_bar(stat = "identity") + facet_wrap(~variable,ncol = 1) + coord_flip() + ylab("Number of metabolites")+xlab("Chemical databases")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                                                                                                                                   panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=14),legend.position = "none",panel.background = element_blank())+scale_fill_npg()

    print(pgS9)
  #---------------------------------Comparing mapping performance rates Figure 3------------------------

    F3<-read.xlsx("CTS_comparisons/CTS-results.xlsx",sheetName = "Mapping-Performance-correct")
  F3$id <-factor(F3$id, levels =  c("MetaFetcheR","MetaboAnalyst","CTS","MS_Targeted"))
    #levels(F3$id)= c("MetaFetcheR","MetaboAnalyst","CTS","MS_Targeted")

    ggplot(F3, aes(x=case,y=value,fill=id)) +
     # geom_bar(stat="identity")  +
     # geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count")
      geom_bar(position="dodge",stat="identity")+ylab("Mapping rate")+xlab("case studies")+guides(fill=guide_legend(title=""))+theme(panel.grid.major = element_blank(),
                                                                                           panel.grid.minor = element_blank(),axis.line = element_line(colour = "grey"),text = element_text(size=20),panel.background = element_blank())+scale_fill_npg()+
      facet_grid(~group)
