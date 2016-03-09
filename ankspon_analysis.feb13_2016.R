#Analysis for Ank Spon data
#based on the biom file 2015-04-24
# JB Feb13, 2016
#################################################################
#required libraries
require(biom)
require(ALDEx2)
require(vegan)
require(ape)
require(ggplot2)
require(GUniFrac)
require(grid)
require(gplots)


#################################################################
#                     Functions Below                           #
#                                                               #
#################################################################
Subsample.Table<-function(OTUTABLE,DEPTH){
  print(paste("Subsampling OTU table to", DEPTH))
  subsampled.OTUTABLE<-t(vegan::rrarefy(t(OTUTABLE),DEPTH)) #expecting the transpose for otu table layout so transpose then transpose back
  subsampled.OTUTABLE<-subsampled.OTUTABLE[rowSums(subsampled.OTUTABLE)>0,] #remove now zero count OTUs
  print(paste("...sampled to",sum(subsampled.OTUTABLE)/ncol(subsampled.OTUTABLE), "with", nrow(subsampled.OTUTABLE), "taxa"))
  return(subsampled.OTUTABLE)
} #subsample with replacement as per vegan implementation, remove otus with under 1 count

Create.Taxa.Lookup<-function(OTUTABLE,TAXONOMY){
  TAX.TABLE<-TAXONOMY
  TAX.TABLE$taxonomy<-rownames(OTUTABLE)

  TAX.TABLE<-apply(TAX.TABLE[,1:ncol(TAX.TABLE)], 2, function(x)(gsub("[^[:alnum:][:space:]_]", "", x))) # remove anything that isn't alpha numeric, a space or an underscore, Green Genes can often have these
  
  #replace the NAs with the appropriate fill in
  TAX.TABLE[is.na(TAX.TABLE[,1]),1]<-"k__"
  TAX.TABLE[is.na(TAX.TABLE[,2]),2]<-"p__"
  TAX.TABLE[is.na(TAX.TABLE[,3]),3]<-"c__"
  TAX.TABLE[is.na(TAX.TABLE[,4]),4]<-"o__"
  TAX.TABLE[is.na(TAX.TABLE[,5]),5]<-"f__"
  TAX.TABLE[is.na(TAX.TABLE[,6]),6]<-"g__"
  TAX.TABLE[is.na(TAX.TABLE[,7]),7]<-"s__"
  TAX.TABLE<-as.data.frame(TAX.TABLE)
  TAX.TABLE$mergedtaxonomy<-paste(TAX.TABLE$taxonomy1,TAX.TABLE$taxonomy2,TAX.TABLE$taxonomy3,TAX.TABLE$taxonomy4,TAX.TABLE$taxonomy5,TAX.TABLE$taxonomy6,TAX.TABLE$taxonomy7,sep=";") #unnecessarily long but crashed when referencing range of columns
  row.names(TAX.TABLE)<-TAX.TABLE$taxonomy
  TAX.TABLE$taxonomy<-NULL
  return(TAX.TABLE)
} #Returns a lookup table of taxonomy where rowname is OTU#, first seven columns are taxonomic levels and 8th column is merged taxonomy. Also fixing the missing values and setting NAs to g__ s__ for example such that every OTUs has these for every rank

Summarize.Taxa<-function(OTUTABLE, TAXONOMY){
  
  SUMMARIZED.TAXA<-list()
  #TaxonomicLevel is table where rowname is OTU, and first column is the merged taxonomy at the given taxonomic level
  for (i in 1:7){
    print(paste("Now summarizing taxonomic level", i, "....."))
    if(i==1){TaxonomicLevel<-cbind(rownames(TAXONOMY),TAXONOMY$taxonomy1)}
    if(i==2){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$taxonomy1,TAXONOMY$taxonomy2,sep=";"))}
    if(i==3){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$taxonomy1,TAXONOMY$taxonomy2,TAXONOMY$taxonomy3,sep=";"))}
    if(i==4){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$taxonomy1,TAXONOMY$taxonomy2,TAXONOMY$taxonomy3,TAXONOMY$taxonomy4,sep=";"))}
    if(i==5){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$taxonomy1,TAXONOMY$taxonomy2,TAXONOMY$taxonomy3,TAXONOMY$taxonomy4,TAXONOMY$taxonomy5,sep=";"))}
    if(i==6){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$taxonomy1,TAXONOMY$taxonomy2,TAXONOMY$taxonomy3,TAXONOMY$taxonomy4,TAXONOMY$taxonomy5,TAXONOMY$taxonomy6,sep=";"))}
    if(i==7){TaxonomicLevel<-cbind(rownames(TAXONOMY),paste(TAXONOMY$taxonomy1,TAXONOMY$taxonomy2,TAXONOMY$taxonomy3,TAXONOMY$taxonomy4,TAXONOMY$taxonomy5,TAXONOMY$taxonomy6,TAXONOMY$taxonomy7,sep=";"))}
    rownames(TaxonomicLevel)<-TaxonomicLevel[,1]
    colnames(TaxonomicLevel)<-c("temp","MergedTaxonomy")

    merged.OTUTABLE<-cbind(as.data.frame(TaxonomicLevel),as.data.frame(OTUTABLE))
    merged.OTUTABLE$temp<-NULL
    Summarized.OTUTABLE<-aggregate(. ~ MergedTaxonomy, data=merged.OTUTABLE, FUN=sum) #aggregate on merged taxonomy ID and sum all the rows collapsed
    rownames(Summarized.OTUTABLE)<-Summarized.OTUTABLE$MergedTaxonomy
    Summarized.OTUTABLE$MergedTaxonomy<-NULL
    SUMMARIZED.TAXA[[length(SUMMARIZED.TAXA)+1]] <- Summarized.OTUTABLE
  }
  taxonomy.OTUTABLE<-OTUTABLE #with taxonomy merged to OTUIDS and all items
  rownames(taxonomy.OTUTABLE)<-paste(row.names(taxonomy.OTUTABLE), TAXONOMY$mergedtaxonomy,sep="|")
  SUMMARIZED.TAXA[[length(SUMMARIZED.TAXA)+1]]<-taxonomy.OTUTABLE
  return(SUMMARIZED.TAXA)
} #Takes an OTU table and a master data frame of OTUs to taxonomy, returns a list of dataframes where the list index is the taxonomic level, the 8th level is OTU

#################################################################
#                     Settings below                            #
#                                                               #
#################################################################
setwd("~/Documents/Research/as_plaques/feb2016/")
metadata<-read.table("as_metadata.txt", header=T, row.names=1, comment.char="", sep="\t")
master.otutable<-read_biom("td_OTU_RDPlineage_fixed.biom") #an unfiltered table (in biom format) exactly as it comes from pipeline of GB Gloor (https://github.com/ggloor/miseq_bin.git)
set.seed(182) #set consistent randomization seed to allow for replication
ascol<-rgb(1,0,0,0.5,maxColorValue = 1)#color to plot AxSpa as
concol<-rgb(0,0,1,0.5,maxColorValue = 1) #color to plot Controls as

#################################################################
#                      Start analysis below                     #
#                                                               #
#################################################################

taxonomy<-observation_metadata(master.otutable) #remove taxonomy ala biom package
master.otutable<-as(biom_data(master.otutable),"matrix") #as per vignette to make workable as matrix
df.taxonomy<-Create.Taxa.Lookup(master.otutable,taxonomy)#a data frame of taxonomies to use as a look up table
summarized.master.otutable<-Summarize.Taxa(master.otutable,df.taxonomy)#make summarized otu tables, returns a list of dataframes. To access: summarized.master.otutable[[taxonomiclevel]], 1-7 correspond to qiimesque L#s, 8th is OTU level


rel.summarized.master.otutable<-list()
for(i in 1:8){
  rel.summarized.master.otutable[[i]]<-apply(summarized.master.otutable[[i]],2, function(x){100*(x/sum(x))})
}#make version of summarized tables in relative abundances (%)
#printed out the tables to double check against qiime's summarize.taxa.py....success!

clr.summarized.master.otutable<-list()
for(p in 1:8){
  #add simpliest prior for log transformation
  clr.summarized.master.otutable[[p]]<-summarized.master.otutable[[p]]+0.5
  clr.summarized.master.otutable[[p]]<-apply(clr.summarized.master.otutable[[p]],2, function(x){log2(x) - mean(log2(x))})
}#make a version of hte table using clr noramlized relative abundances

########Do Alpha diversity analysis

sub.OTUlevel<-Subsample.Table(master.otutable,min(colSums(master.otutable))) #subsample to lowest number of reads which is 9064

#Do alpha diversity metrics with vegan
shannon<-vegan::diversity(t(sub.OTUlevel), index ="shannon") #from vegan
  metadata_alpha<-merge(metadata,as.data.frame(shannon),by="row.names",all=T)
  rownames(metadata_alpha)<-metadata_alpha$Row.names
  metadata_alpha$Row.names<-NULL
simpson<-vegan::diversity(t(sub.OTUlevel), index ="simpson") #from vegan
  metadata_alpha<-merge(metadata_alpha,as.data.frame(simpson),by="row.names",all=T)
  rownames(metadata_alpha)<-metadata_alpha$Row.names
  metadata_alpha$Row.names<-NULL
chao<-vegan::estimateR(t(sub.OTUlevel)) #get chao1 and ace from vegan
  metadata_alpha<-merge(metadata_alpha,as.data.frame(t(chao)),by="row.names",all=T)

rownames(metadata_alpha)<-metadata_alpha$Row.names
metadata_alpha$Row.names<-NULL

metadata_alpha$group<-gsub("as","Ankylosing Spondylitis", metadata_alpha$group) #recode variables to be more informative
metadata_alpha$group<-gsub("control","Healthy Control", metadata_alpha$group)  

shannonplot<-ggplot(metadata_alpha, aes(x=group, y=shannon, color=group)) + geom_boxplot()+geom_jitter(position=position_jitter(0.2)) +labs(x=NULL, y="Shannon's Diversity Index (H)") +theme_bw() + scale_color_manual(values=c(ascol,concol))  + theme(legend.position = "none", axis.text.x = element_blank())
simpsonplot<-ggplot(metadata_alpha, aes(x=group, y=simpson, color=group)) + geom_boxplot()+geom_jitter(position=position_jitter(0.2)) +labs(x=NULL, y="Simpsons's Diversity Index (D)") +theme_bw() + scale_color_manual(values=c(ascol,concol))  + theme(legend.position = "none", axis.text.x = element_blank())
chaoplot<-ggplot(metadata_alpha, aes(x=group, y=S.chao1, color=group)) + geom_boxplot()+geom_jitter(position=position_jitter(0.2)) +labs(x=NULL, y="Chao1 Richness Estimate") +theme_bw() + scale_color_manual(values=c(ascol,concol))  + theme(legend.position = "none", axis.text.x = element_blank())

wilcox.test(subset(metadata_alpha, group=="Ankylosing Spondylitis")$shannon, subset(metadata_alpha, group=="Healthy Control")$shannon)
wilcox.test(subset(metadata_alpha, group=="Ankylosing Spondylitis")$simpson, subset(metadata_alpha, group=="Healthy Control")$simpson)
wilcox.test(subset(metadata_alpha, group=="Ankylosing Spondylitis")$S.chao1, subset(metadata_alpha, group=="Healthy Control")$S.chao1)

########Do Beta diversity analysis

tree<-ape::read.tree("as_otu.tree") #result of Muscle alignment and make_phylogeny.py with fasttree and midpoint rooting from MACQIIME 1.9.0-20140227
unifracs <- GUniFrac::GUniFrac(t(sub.OTUlevel),tree, alpha=c(0, 0.5, 1))$unifracs

weighted.UniFrac <- unifracs[, , "d_1"] # Weighted UniFrac as per ?GUniFrac
unweighted.UniFrac <- unifracs[, , "d_UW"] # Unweighted UniFrac
braycurtis<-vegan::vegdist(t(sub.OTUlevel), method="bray", binary=FALSE, diag=FALSE, upper=FALSE) 
braycurtis<-as.matrix(braycurtis)

vegan::anosim(weighted.UniFrac, metadata[rownames(weighted.UniFrac),]$group, permutations=999)
vegan::anosim(unweighted.UniFrac, metadata[rownames(unweighted.UniFrac),]$group, permutations=999)
vegan::anosim(braycurtis, metadata[rownames(braycurtis),]$group, permutations=999)

metadata<-metadata[rownames(weighted.UniFrac),] #filter metadata to only include samples in analysis and reorder

pco.braycurtis<-ape::pcoa(braycurtis) #using the principal coordinants analysis of the ape package
pco.weighted.UniFrac<-ape::pcoa(weighted.UniFrac)
pco.unweighted.UniFrac<-ape::pcoa(unweighted.UniFrac)

plot.pco.braycurtis<-merge(as.data.frame(pco.braycurtis$vectors), metadata[rownames(pco.braycurtis$vectors),], by="row.names", all.x=T)
plot.pco.weighted.UniFrac<-merge(as.data.frame(pco.weighted.UniFrac$vectors), metadata[rownames(pco.weighted.UniFrac$vectors),], by="row.names", all.x=T)
plot.pco.unweighted.UniFrac<-merge(as.data.frame(pco.unweighted.UniFrac$vectors), metadata[rownames(pco.unweighted.UniFrac$vectors),], by="row.names", all.x=T)


pc.braycurtis<-ggplot(plot.pco.braycurtis, aes(x=Axis.1, y=Axis.2, color=group)) + geom_point()+ scale_colour_manual('Group', values = c(ascol,concol)) + theme_bw() + ggtitle("Bray Curtis") + theme(legend.position = "none") + coord_fixed(ratio=1) + xlab(paste("PCo1:", round((100*pco.braycurtis$values$Eigenvalues/sum(pco.braycurtis$values$Eigenvalues))[1],1), "% Variation Explained")) + ylab(paste("PCo2:", round((100*pco.braycurtis$values$Eigenvalues/sum(pco.braycurtis$values$Eigenvalues))[2],1), "% Variation Explained"))
pc.weighted.UniFrac<-ggplot(plot.pco.weighted.UniFrac, aes(x=Axis.1, y=Axis.2, color=group)) + geom_point()+ scale_colour_manual('Group', values = c(ascol,concol)) + theme_bw() + ggtitle("Weighted UniFrac") + theme(legend.position = "none") + coord_fixed(ratio=1) + xlab(paste("PCo1:", round((100*pco.weighted.UniFrac$values$Eigenvalues/sum(pco.weighted.UniFrac$values$Eigenvalues))[1],1), "% Variation Explained")) + ylab(paste("PCo2:", round((100*pco.weighted.UniFrac$values$Eigenvalues/sum(pco.weighted.UniFrac$values$Eigenvalues))[2],1), "% Variation Explained"))
pc.unweighted.UniFrac<-ggplot(plot.pco.unweighted.UniFrac, aes(x=Axis.1, y=Axis.2, color=group)) + geom_point()+ scale_colour_manual('Group', values = c(ascol,concol)) + theme_bw() + ggtitle("Unweighted UniFrac") + theme(legend.position = "none") + coord_fixed(ratio=1) + xlab(paste("PCo1:", round((100*pco.unweighted.UniFrac$values$Eigenvalues/sum(pco.unweighted.UniFrac$values$Eigenvalues))[1],1), "% Variation Explained")) + ylab(paste("PCo2:", round((100*pco.unweighted.UniFrac$values$Eigenvalues/sum(pco.unweighted.UniFrac$values$Eigenvalues))[2],1), "% Variation Explained"))

####now plot
pdf("Diversity_analysis.pdf",height=5,width=8,useDingbats=FALSE)
  pushViewport(viewport(layout = grid.layout(2, 3)))
  print(shannonplot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(simpsonplot, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(chaoplot, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
  print(pc.braycurtis, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(pc.weighted.UniFrac, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  print(pc.unweighted.UniFrac, vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
dev.off()

##################################################################################
#Make a new Heatmap at the family level, alter i to change taxonomic level plotted. Differs from publication in aesthetics, log transformation, and row ordering

as.samps<-rownames(subset(metadata, group=="as"))
con.samps<-rownames(subset(metadata, group=="control"))
i=5 #plot family level, can be modified to plot any level with 8 OTU level

#UPGMA clustering both on euclidian distances for plotting purposes, similar to default parameters of heatmap.2
asclust<-colnames(rel.summarized.master.otutable[[i]][,as.samps])[hclust(dist(t(rel.summarized.master.otutable[[i]][,as.samps])))$order]
conclust<-colnames(rel.summarized.master.otutable[[i]][,con.samps])[hclust(dist(t(rel.summarized.master.otutable[[i]][,con.samps])))$order]

toplot<-cbind(rel.summarized.master.otutable[[i]][,asclust], rep(NA,nrow(rel.summarized.master.otutable[[i]])), rel.summarized.master.otutable[[i]][,conclust])
toplot=toplot+0.01
toplot<-log10(toplot)
pdf("as_genus_heatmap.pdf", height=6,width=8, useDingbats = F)
heatmap.2(toplot,
          dendrogram="row",
          Rowv=FALSE,
          Colv=FALSE,
          revC=FALSE,
          trace="none",
          keysize="1",
          density.info="none",
          sepwidth=c(0.05,0.05),
          col=colorRampPalette(c("black","blue","cyan","green","yellow","red"))(100),
          margins=c(1,30),
          labRow=gsub("Bacteria;","",row.names(toplot)))
dev.off()


##########do aldex analysis for cross section between AnkSpon and Control samples
as.METADATA<-subset(metadata,group=="as")
con.METADATA<-subset(metadata,group=="control")

OTUTABLE<-summarized.master.otutable
complete.aldex<-list()  
  for(i in 2:length(OTUTABLE)){
    as.OTUTABLE<-OTUTABLE[[i]][,rownames(as.METADATA)]
    con.OTUTABLE<-OTUTABLE[[i]][,rownames(con.METADATA)]
    
    foraldex<-merge(as.OTUTABLE,con.OTUTABLE, by="row.names", all=TRUE)
    row.names(foraldex)<-foraldex$Row.names
    foraldex$Row.names<-NULL
    
    conds<-c(rep("AnkSpon",ncol(as.OTUTABLE)),rep("HealthyControl",ncol(con.OTUTABLE)))
    
    aldex.clr<-aldex.clr(foraldex, mc.samples=256)
    otutable.aldex<-aldex.ttest(aldex.clr, conds, paired.test=FALSE)
    otutable.aldex.effect<-aldex.effect(aldex.clr, conds, useMC=TRUE, include.sample.summary = TRUE)
    complete.aldex[[i]]<-merge(otutable.aldex,otutable.aldex.effect, by="row.names", all=TRUE)
    rownames(complete.aldex[[i]])<-complete.aldex[[i]]$Row.names
    complete.aldex[[i]]$Row.names<-NULL
  }
    pdf("aldexplot.pdf",height=4,width=6,useDingbats = F)
    plot(complete.aldex[[2]]$diff.win, complete.aldex$diff.btw, xlab="Median Log2 Within Condition Difference", ylab="Median Log2 Between Condition Difference", col=rgb(0,0,0,0.6,maxColorValue=1), pch=16, cex=0.5,xlim=c(0,10),ylim=c(-10,10))
    abline(0,1, col="grey", lty=2, lwd=1)
    abline(0,-1, col="grey", lty=2, lwd=1)
    
    for (j in 3:length(complete.aldex)){
      points(complete.aldex[[j]]$diff.win,complete.aldex[[j]]$diff.btw, col=rgb(0,0,0,0.6,maxColorValue=1), pch=16, cex=0.5)
    }
    dev.off()
    
    for (p in 2:length(complete.aldex)){
    write.table(complete.aldex[[p]],paste("L",p,"aldex.table.txt",sep=""), quote=F, sep='\t', col.names=NA)
    }
    
    #collect all comparisons and readjust pvalues for all comparisons, individual outputs are for single comparison
    finalp<-rbind(complete.aldex[[2]],complete.aldex[[3]],complete.aldex[[4]],complete.aldex[[5]],complete.aldex[[6]],complete.aldex[[7]],complete.aldex[[8]])
    finalp$we.eBH<-p.adjust(finalp$we.ep,method="BH", n=length(finalp$we.ep))
    finalp$wi.eBH <-p.adjust(finalp$wi.ep,method="BH", n=length(finalp$wi.ep))
    
    write.table(finalp,paste("all_comparisons_corrected.aldex.table.txt",sep=""), quote=F, sep='\t', col.names=NA)
    message("There are: ", sum(finalp$we.eBH<0.1), " significant features with a corrected-p<0.1!")

  #Now do correlations with metadata using CLR-normalized data
    meta.to.compare<-c("npd4plus","BASDAI","ASQoL","CRP") #list of metadata to use for correlation
    corr.results<-data.frame(c("OTU",rownames(clr.summarized.master.otutable[[8]])))
    for(y in 1:length(meta.to.compare)){
      pvals<-paste(meta.to.compare[y],".p.value",sep="")
      rvals<-paste(meta.to.compare[y],".r.value",sep="")
      #rownames(metadata)==colnames(clr.summarized.master.otutable[[8]]) #double checking correct ordering of sampleIDs
      for (x in 1:nrow(clr.summarized.master.otutable[[8]])){
       corr<-cor.test(clr.summarized.master.otutable[[8]][x,],metadata[,meta.to.compare[y]], method="pearson")
        r<-cor(clr.summarized.master.otutable[[8]][x,],metadata[,meta.to.compare[y]], method="pearson", use="complete.obs")
        pvals<-append(pvals,corr$p.value)
        rvals<-append(rvals,r)
      }
      pvals<-c(pvals[1],p.adjust(pvals[2:length(pvals)], method="BH", n=length(meta.to.compare)*(length(pvals)-1))) #adjust for all comparisons made
      #pvals<-c(pvals[1],p.adjust(pvals[2:length(pvals)], method="BH", n=length(pvals)-1)) #adjust within a single regression group
      corr.results<-cbind(corr.results,pvals,rvals)
    }
    
    colnames(corr.results)<-t(corr.results)[,1]
    corr.results<-corr.results[-1,]
    
    write.table(corr.results,"correlationanalysis.txt", sep='\t', col.names=NA, quote=F)
    
    
    t<-read.table("correlationanalysis.txt", sep='\t', header=T, row.names=1)
    message("There are: ", sum(t$npd4plus.p.value<0.1), " significant features with a corrected-p<0.1 for npd4plus!")
    message("There are: ", sum(t$BASDAI.p.value<0.1), " significant features with a corrected-p<0.1 for BASDAI!")
    message("There are: ", sum(t$CRP.p.value<0.1), " significant features with a corrected-p<0.1 for CRP!")
    message("There are: ", sum(t$ASQoL.p.value<0.1), " significant features with a corrected-p<0.1 for ASQoL!")
    
    message("Analysis complete")
    
    