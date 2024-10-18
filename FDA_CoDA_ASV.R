#FDA project
#Analysis of microbiomes and mycobiomes using CoDa approach at ASV level for two years
#Laura Rolon & Marysabel Mendez-Acevedo
#Last updated: 10/06/2023 MLR

#Load packages
library(ape)
library(vegan)
library(ggplot2)
library(phyloseq)
library(tidyr)
library(dplyr)
library(compositions)
library(zCompositions)
library(viridis)
library(readxl)
library(pairwiseAdonis)
library(psych)
library(svglite)
library(decontam)

set.seed(336)

#Set working directory to where files are located
setwd("G:/My Drive/Penn State/Research/File for R/FDA_PDA/TwoYears")

#### IMPORT DATA ####

#Import data - 16s 
asvs_16s<-read.csv('ASV_16s.csv', header = TRUE, row.names = 1)
taxon_16s<-read.csv('Taxon_16s.csv', header = TRUE, row.names = 1)
metadata_16s<-as.data.frame(read_excel('Results_FDAPDA.xlsx', sheet="Metadata_16S", col_names = TRUE))
row.names(metadata_16s)<-metadata_16s$SampleID


#Import ASV table - ITS
asvs_ITS<-read.csv('ASV_ITS.csv', header = TRUE, row.names = 1)
taxon_ITS<-read.csv('Taxon_ITS.csv', header = TRUE, row.names = 1)
metadata_ITS<-as.data.frame(read_excel('Results_FDAPDA.xlsx', sheet="Metadata_ITS", col_names = TRUE))
row.names(metadata_ITS)<-metadata_ITS$SampleID

#Clean up Taxon tables
#Add '_unclassified' marker to NAs in taxon table
taxon_16s$Phylum<-ifelse(is.na(taxon_16s$Phylum), paste(taxon_16s$Kingdom, "unclassified", sep = '_'), taxon_16s$Phylum)
taxon_16s$Class<-ifelse(is.na(taxon_16s$Class), paste(taxon_16s$Phylum, "unclassified", sep = '_'), taxon_16s$Class)
taxon_16s$Order<-ifelse(is.na(taxon_16s$Order), paste(taxon_16s$Class, "unclassified", sep = '_'), taxon_16s$Order)
taxon_16s$Family<-ifelse(is.na(taxon_16s$Family), paste(taxon_16s$Order, "unclassified", sep = '_'), taxon_16s$Family)
taxon_16s$Genus<-ifelse(is.na(taxon_16s$Genus), paste(taxon_16s$Family, "unclassified", sep = '_'), taxon_16s$Genus)
taxon_16s$Species<-ifelse(is.na(taxon_16s$Species), paste(taxon_16s$Genus, "unclassified", sep = '_'), taxon_16s$Species)

taxon_ITS$Phylum<-ifelse(is.na(taxon_ITS$Phylum), paste(taxon_ITS$Kingdom, "unclassified", sep = '_'), taxon_ITS$Phylum)
taxon_ITS$Class<-ifelse(is.na(taxon_ITS$Class), paste(taxon_ITS$Phylum, "unclassified", sep = '_'), taxon_ITS$Class)
taxon_ITS$Order<-ifelse(is.na(taxon_ITS$Order), paste(taxon_ITS$Class, "unclassified", sep = '_'), taxon_ITS$Order)
taxon_ITS$Family<-ifelse(is.na(taxon_ITS$Family), paste(taxon_ITS$Order, "unclassified", sep = '_'), taxon_ITS$Family)
taxon_ITS$Genus<-ifelse(is.na(taxon_ITS$Genus), paste(taxon_ITS$Family, "unclassified", sep = '_'), taxon_ITS$Genus)
taxon_ITS$Species<-ifelse(is.na(taxon_ITS$Species), paste(taxon_ITS$Genus, "unclassified", sep = '_'), taxon_ITS$Species)

#Remove extra _unclassified
taxon_16s$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Class)
taxon_16s$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Order)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Family)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Genus)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_16s$Species)
taxon_16s$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_16s$Species)

taxon_ITS$Class<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Class)
taxon_ITS$Order<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Order)
taxon_ITS$Order<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Order)
taxon_ITS$Family<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Family)
taxon_ITS$Family<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Family)
taxon_ITS$Family<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Family)
taxon_ITS$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Genus)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified_unclassified", '_unclassified', taxon_ITS$Species)
taxon_ITS$Species<-gsub("_unclassified_unclassified", '_unclassified', taxon_ITS$Species)

#Remove taxonomic rank from ITS taxon table
taxon_ITS$Kingdom<-gsub("k__","", taxon_ITS$Kingdom)
taxon_ITS$Phylum<-gsub("k__","", taxon_ITS$Phylum)
taxon_ITS$Phylum<-gsub("p__","", taxon_ITS$Phylum)
taxon_ITS$Class<-gsub("k__","", taxon_ITS$Class)
taxon_ITS$Class<-gsub("p__","", taxon_ITS$Class)
taxon_ITS$Class<-gsub("c__","", taxon_ITS$Class)
taxon_ITS$Order<-gsub("k__","", taxon_ITS$Order)
taxon_ITS$Order<-gsub("p__","", taxon_ITS$Order)
taxon_ITS$Order<-gsub("c__","", taxon_ITS$Order)
taxon_ITS$Order<-gsub("o__","", taxon_ITS$Order)
taxon_ITS$Family<-gsub("k__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("p__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("c__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("o__","", taxon_ITS$Family)
taxon_ITS$Family<-gsub("f__","", taxon_ITS$Family)
taxon_ITS$Genus<-gsub("k__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("p__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("c__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("o__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("f__","", taxon_ITS$Genus)
taxon_ITS$Genus<-gsub("g__","", taxon_ITS$Genus)
taxon_ITS$Species<-gsub("k__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("p__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("c__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("o__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("f__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("g__","", taxon_ITS$Species)
taxon_ITS$Species<-gsub("s__","", taxon_ITS$Species)

#Convert asv and taxon tables to matrix
asvs_16s<-as.matrix(asvs_16s)
taxon_16s<-as.matrix(taxon_16s)

asvs_ITS<-as.matrix(asvs_ITS)
taxon_ITS<-as.matrix(taxon_ITS)

#Make phyloseq object
ps_16s<-phyloseq(otu_table(asvs_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s))
ps_ITS<-phyloseq(otu_table(asvs_ITS, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata_ITS))

#Remove Chloroplast and Mitochondria reads from ASV table
physeq_16s <- ps_16s %>%  subset_taxa( Order!="Chloroplast" | is.na(Order) )
physeq_16s <- physeq_16s %>% subset_taxa( Family!= "Mitochondria" | is.na("Family"))

#Get clean ASV table from phyloseq object
asv_16s<-as.data.frame(t(otu_table(physeq_16s)))
tail(rowSums(asv_16s))

asv_ITS<-as.data.frame(t(otu_table(ps_ITS)))
tail(rowSums(asv_ITS))

#Get Taxon table from phyloseq object
taxon.16s<-as.matrix(tax_table(ps_16s))
taxon.ITS<-as.matrix(tax_table(ps_ITS))

#Remove ASVs with zero counts in all samples
asv_16s<-asv_16s[ which(rowSums(asv_16s)>0),]
asv_16s<-t(asv_16s)

asv_ITS<-asv_ITS[ which(rowSums(asv_ITS)>0),]
asv_ITS<-t(asv_ITS)


####  Identify ASVs that match positive control strains ####
#Need to do this step before decontaminating data because most of the ASVs in the positive controls are not found in dataset and decontam will consider them as contaminants

#Obtain positive control declared composition and sequences from Zymo website: https://s3.amazonaws.com/zymofiles/BioPool/ZymoBIOMICS.STD.refseq.v2.zip.  
#Match all ASVs of Pseudomonas, E. coli, Salmonella, Lactobacillus, Enterococcus, Staphylococcus, Listeria and Bacillus to Zymo reference genomes using Mega11 to identify ASVs that are organisms in positive control
#Select ASVs that match with 0 SNPs to reference seqs from Zymo

## 16S ##
#ASV31 - Bacillus
#ASV62 - Enterococcus
#ASV47 - Escherichia
#ASV36 - Listeria
#ASV90 - Pseudomonas
#ASV61/226 Salmonella
#ASV37 - Staphylococcus
#ASV66 - Limosilactobacillus / formerly Lactobacillus

## ITS ## 
#Cannot align because Zymo only provides 18S and mitochondrial sequences

#Make Phyloseq with cleaned ASV
phyloseq_16s<-phyloseq(otu_table(asv_16s, taxa_are_rows = FALSE), tax_table(taxon.16s), sample_data(metadata_16s))

phyloseq_PC<-subset_samples(phyloseq_16s, Control =="Positive")
phyloseq_PCDNA<-subset_samples(phyloseq_16s, Control =="PositiveDNA")

phyloseq_PC_clean<-subset_taxa(phyloseq_PC, rownames(tax_table(phyloseq_PC))=="ASV31"|
                                 rownames(tax_table(phyloseq_PC))=="ASV62"|
                                 rownames(tax_table(phyloseq_PC))=="ASV47"|
                                 rownames(tax_table(phyloseq_PC))=="ASV36"|
                                 rownames(tax_table(phyloseq_PC))=="ASV90"|
                                 rownames(tax_table(phyloseq_PC))=="ASV61"|
                                 rownames(tax_table(phyloseq_PC))=="ASV37"|
                                 rownames(tax_table(phyloseq_PC))=="ASV66")
phyloseq_PCDNA_clean<-subset_taxa(phyloseq_PCDNA, rownames(tax_table(phyloseq_PCDNA))=="ASV31"|
                                                                       rownames(tax_table(phyloseq_PCDNA))=="ASV62"|
                                                                       rownames(tax_table(phyloseq_PCDNA))=="ASV47"|
                                                                       rownames(tax_table(phyloseq_PCDNA))=="ASV36"|
                                                                       rownames(tax_table(phyloseq_PCDNA))=="ASV90"|
                                                                       rownames(tax_table(phyloseq_PCDNA))=="ASV61"|
                                                                       rownames(tax_table(phyloseq_PCDNA))=="ASV37"|
                                                                       rownames(tax_table(phyloseq_PCDNA))=="ASV66")

#See positive extraction control by rep
PC_RA_rep<-phyloseq_PC_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()

barplot_PC_reps<-ggplot(PC_RA_rep, aes(x=SampleID,y=Abundance, fill=Genus))+
  geom_bar(stat='identity', color='black')+facet_grid(.~Year, scales = "free_x")+
  geom_text(aes(label=round(Abundance, digits=1)), position = position_stack(vjust=0.5), size=4)+ylab("Relative abundance (%)")+
  ggtitle("PC barplot by rep")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))

PCDNA_RA<-phyloseq_PCDNA_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()

PC_Zymo<-data.frame(OTU=c("ASV31","ASV62","ASV47","ASV36","ASV90","ASV61","ASV37","ASV66"),
                    Sample=rep("Zymo",8),
                    Genus=c("Bacillus","Enterococcus","Escherichia-Shigella","Listeria","Pseudomonas","Salmonella","Staphylococcus","Limosilactobacillus"),
                    Abundance=c(17.4,9.9,10.1,14.1,4.2,10.4,15.5,18.4))

PC_data_all<-bind_rows(PC_RA_rep,PCDNA_RA,PC_Zymo)
PC_data_all$Abundance<-round(PC_data_all$Abundance, digits = 1)
PC_data_all<-PC_data_all[order(PC_data_all$Sample),]
PC_data_all$SampleOrder<-c(rep(11,8),rep(5,8),rep(2,8),rep(8,8),rep(6,8),rep(3,8),rep(9,8),rep(7,8),rep(4,8),rep(10,8),rep(1,8))

barplot_PC_reps<-ggplot(PC_data_all, aes(x=reorder(Sample,SampleOrder),y=Abundance, fill=Genus))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=Abundance), position = position_stack(vjust=0.5), size=4)+  ylab("Relative abundance (%)")+
  ggtitle("PC barplot")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Positive controls_reps.png", plot =barplot_PC_reps, device="png", width=20, height=10, units="in",dpi=600)
ggsave("Positive controls_reps.svg", plot =barplot_PC_reps, device="svg", width=20, height=10, units="in",dpi=600)

#Calculate relative abundance and average across replicates
PC_RA<-phyloseq_PC_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()%>%
  group_by(OTU, Genus, Control, Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

PCDNA_RA<-phyloseq_PCDNA_clean%>%
  transform_sample_counts(function(x) x/sum(x)*100)%>%
  psmelt()%>%
  group_by(OTU, Genus, Control,Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

PC_Zymo<-data.frame(OTU=c("ASV31","ASV62","ASV47","ASV36","ASV90","ASV61","ASV37","ASV66"),
                    Control=rep("Zymo",8),
                    Genus=c("Bacillus","Enterococcus","Escherichia-Shigella","Listeria","Pseudomonas","Salmonella","Staphylococcus","Limosilactobacillus"),
                    Mean=c(17.4,9.9,10.1,14.1,4.2,10.4,15.5,18.4))

PC_data<-bind_rows(PC_RA,PCDNA_RA,PC_Zymo)
PC_data$Mean<-round(PC_data$Mean, digits = 1)
PC_data<-PC_data[order(PC_data$Control),]
PC_data$Control_Year<-interaction(PC_data$Control,PC_data$Year)
PC_data$Order<-c(rep(c(2,3),8),rep(c(4,5),8),rep(1,8))

#Plot with average 
barplot_PC_all<-ggplot(PC_data, aes(x=reorder(Control_Year,Order),y=Mean, fill=Genus))+
  geom_bar(stat='identity', color='black')+
  geom_text(aes(label=Mean), position = position_stack(vjust=0.5), size=4)+  ylab("Relative abundance (%)")+
  ggtitle("PC barplot")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'),
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Positive controls.png", plot =barplot_PC_all, device="png", width=20, height=10, units="in",dpi=600)
ggsave("Positive controls.svg", plot =barplot_PC_all, device="svg", width=20, height=10, units="in",dpi=600)



#### Check negative control and decontaminate ASVs ####
phyloseq_16s<-phyloseq(otu_table(asv_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s))
phyloseq_ITS<-phyloseq(otu_table(asv_ITS, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata_ITS))

#Subset negative controls
phyloseq_NC<-subset_samples(phyloseq_16s, Control =="Negative")
phyloseq_NC_ITS<-subset_samples(phyloseq_ITS, Control =="Negative")

#Make long table
phyloseq.NC<-psmelt(phyloseq_NC)
phyloseq.NC_ITS<-psmelt(phyloseq_NC_ITS)

#Remove rows with less than 1000 reads
phyloseq.NC<-subset(phyloseq.NC, Abundance>1000)
phyloseq.NC_ITS<-subset(phyloseq.NC_ITS, Abundance>1000)

#Plot reads of control by  NC
barplot_NC<-ggplot(phyloseq.NC, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=SampleID))+
  geom_bar(stat='identity', color='black')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Negative control read counts")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_NC.png", plot =barplot_NC, device="png", width=20, height=10, units="in",dpi=600)

barplot_NC_ITS<-ggplot(phyloseq.NC_ITS, aes(x=reorder(OTU, desc(Abundance)),y=Abundance, fill=SampleID))+
  geom_bar(stat='identity', color='black')+ 
  geom_text(aes(label=Genus, angle=90, hjust=0, size=13))+  ylab("Reads (x10,000)")+
  scale_y_continuous(breaks= seq(0,75000, 5000), labels = function(x){x/10000}, limits=c(0,75000)) + 
  ggtitle("Negative control read counts")+xlab("")+
  theme(axis.text.y=element_text(color='black', size=13), axis.ticks=element_line(color='black'),
        axis.text.x = element_text(angle=90, size=13, color='black'), axis.title.y = element_text(color='black', size=13))+
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 1)) +
  theme(legend.position = "none")+
  theme(panel.background = element_rect(fill="grey99", color =NA),plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
ggsave("Controls_NC_TS.png", plot =barplot_NC_ITS, device="png", width=20, height=10, units="in",dpi=600)

#### Decontaminate sequencing data ####
metadata_16s$Year<-as.factor(metadata_16s$Year)
metadata_ITS$Year<-as.factor(metadata_ITS$Year)

#Detect contaminants with prevalence method (default threshold=0.1)
sample_data(phyloseq_16s)$is.neg <- sample_data(phyloseq_16s)$Control == "Negative"
contamdf.prev_16s <- isContaminant(phyloseq_16s, method="prevalence", neg="is.neg", threshold=0.1, batch="Year")
table(contamdf.prev_16s$contaminant)
which(contamdf.prev_16s$contaminant) 

sample_data(phyloseq_ITS)$is.neg <- sample_data(phyloseq_ITS)$Control == "Negative"
contamdf.prev_ITS <- isContaminant(phyloseq_ITS, method="prevalence", neg="is.neg", threshold=0.1, batch="Year")
table(contamdf.prev_ITS$contaminant)
which(contamdf.prev_ITS$contaminant) 


#Remove identified contaminants (prevalence method) from ASV table
phyloseq.noncontam_16s <- prune_taxa(!contamdf.prev_16s$contaminant, phyloseq_16s)
phyloseq.noncontam_16s

phyloseq.noncontam_ITS <- prune_taxa(!contamdf.prev_ITS$contaminant, phyloseq_ITS)
phyloseq.noncontam_ITS

#Remove negative and positive controls from ASV table
phyloseq_16s_clean<-subset_samples(phyloseq.noncontam_16s, Control=="No" )
phyloseq_ITS_clean<-subset_samples(phyloseq.noncontam_ITS, Control=="No" )

#Get metadata 
metadata_16s_clean<-subset(metadata_16s, Control== "No")
metadata_ITS_clean<-subset(metadata_ITS, Control== "No")

#Get ASV table from phyloseq object
asv_16s_clean<-t(as.data.frame(otu_table(phyloseq_16s_clean)))
tail(rowSums(asv_16s_clean))

asv_ITS_clean<-t(as.data.frame(otu_table(phyloseq_ITS_clean)))
tail(rowSums(asv_ITS_clean))

#Remove ASVs with zero counts in all samples
asv_16s_clean<-asv_16s_clean[ which(rowSums(asv_16s_clean)>0),]
asv_ITS_clean<-asv_ITS_clean[ which(rowSums(asv_ITS_clean)>0),]


#### COMPOSITIONAL ANALYSIS OF MICROBIOME AT ASV LEVEL ####
#Based on Microbiome Analysis in R. Chap 10.

#Step 1: Convert ASV table to appropriate format
#Following step requires samples on rows and ASVs in columns
head(t(asv_16s_clean)) 
head(t(asv_ITS_clean))

#Step 2: Replace zero values before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16s<-t(cmultRepl(t(asv_16s_clean), label=0, method="CZM", output="p-counts")) #No. adjusted imputations:  6047 
asv.n0_ITS<-t(cmultRepl(t(asv_ITS_clean), label=0, method="CZM", output="p-counts")) # No. adjusted imputations:  44712 

#Note: Check the output to make sure there are no negative numbers. If samples or ASV are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv_n0_16s<-ifelse(asv.n0_16s < 0, asv.n0_16s*(-1), asv.n0_16s)
asv_n0_ITS<-ifelse(asv.n0_ITS < 0, asv.n0_ITS*(-1), asv.n0_ITS)

#output table needs to have samples in columns and ASVs in rows
head(asv_n0_16s) 
head(asv_n0_ITS)

#Step 3: Convert data to proportions
asv.n0_16s_prop<-apply(asv_n0_16s, 2, function(x) {x/sum(x)})
asv.n0_ITS_prop<-apply(asv_n0_ITS, 2, function(x) {x/sum(x)})

#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_16s_prop_f<-asv_n0_16s[apply(asv.n0_16s_prop, 1, min) > 0.000001, ]
asv.n0_ITS_prop_f<-asv_n0_ITS[apply(asv.n0_ITS_prop, 1, min) > 0.000001, ]

#Check that samples are on columns and ASVs in rows
head(asv.n0_16s_prop_f) 
head(asv.n0_ITS_prop_f)

#Step 5: perform CLR transformation
asv.n0.clr_16s<-t(apply(asv.n0_16s_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITS<-t(apply(asv.n0_ITS_prop_f, 2, function(x){log(x)-mean(log(x))}))

#Check output table. Samples should be in rows and ASVs in columns
head(asv.n0.clr_16s) 
head(asv.n0.clr_ITS)

#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16s<-prcomp(asv.n0.clr_16s)
pc.clr_ITS<-prcomp(asv.n0.clr_ITS)

png("Screeplot - PCA.png", width = 400, height = 300, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(1,2))
screeplot(pc.clr_16s, type='barplot', main="Bacteria")
screeplot(pc.clr_ITS, type='barplot', main="Fungi ")
dev.off()

#Calculate total variance of the data
mvar.clr_16s<-mvar(asv.n0.clr_16s)
mvar.clr_ITS<-mvar(asv.n0.clr_ITS)

#Display results - 16s 
row_16s<-rownames(asv.n0.clr_16s) #Make vector with sample names
pc_out_16s<-as.data.frame(pc.clr_16s$x[,1:3]) #Get PC1 and PC2 
pc_out_meta_16s<-as.data.frame(bind_cols(pc_out_16s,metadata_16s_clean)) #Add metadata information. Make sure metadata is organized alphabetically.
row.names(pc_out_meta_16s)<-row_16s #Add rownames to dataframe
pc_out_meta_16s$Facility<-as.factor(pc_out_meta_16s$Facility)
pc_out_meta_16s$Time<-as.factor(pc_out_meta_16s$Time)
pc_out_meta_16s$Treatment<-as.factor(pc_out_meta_16s$Treatment)

loadings_16s<-as.data.frame(pc.clr_16s$rotation)

# Make PCA plot - First 2 axis, color by facility and shape by Year
PCA_16s <- ggplot(pc_out_meta_16s, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
    #scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16s$sdev[1]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16s$sdev[2]^2/mvar.clr_16s*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria", subtitle = "PCA by Facility and Treatment")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_16s
ggsave("PCA_Bacteria_ASV.png", plot =PCA_16s, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Bacteria_ASV.svg", plot =PCA_16s, device="svg", width=6, height=5, units="in",dpi=600)




#Display results -ITS
row_ITS<-rownames(asv.n0.clr_ITS) #Make vector with sample names
pc_out_ITS<-as.data.frame(pc.clr_ITS$x[,1:2]) #Get PC1 and PC2
pc_out_meta_ITS<-as.data.frame(bind_cols(pc_out_ITS,metadata_ITS_clean)) #Add metadata information
row.names(pc_out_meta_ITS)<-row_ITS #Add rownames to dataframe
pc_out_meta_ITS$Facility<-as.factor(pc_out_meta_ITS$Facility)
pc_out_meta_ITS$Time<-as.factor(pc_out_meta_ITS$Time)
pc_out_meta_ITS$Treatment<-as.factor(pc_out_meta_ITS$Treatment)

# Make PCA plot - First 2 axis, color by facility and shape by Year
PCA_ITS <- ggplot(pc_out_meta_ITS, aes(x=PC1,y=PC2, color=Facility, shape=Year))+
  geom_point(size=3)+
  #scale_shape_manual(values = c(15,12))+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=13, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITS$sdev[1]^2/mvar.clr_ITS*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITS$sdev[2]^2/mvar.clr_ITS*100, digits=1), "%", sep="")) +
  ggtitle("Fungi", subtitle = "PCA by Facility and Season")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')
PCA_ITS
ggsave("PCA_Fungi_ASV.png", plot =PCA_ITS, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Fungi_ASV.svg", plot =PCA_ITS, device="svg", width=6, height=5, units="in",dpi=600)

#### Try package Microviz ####
library(microViz)

ord_explore(phyloseq_16s_clean) #This will sent to a webpage to set the ordination options. Copy code when satisfied with result

PC_16s_loadings<-phyloseq_16s_clean %>%
  tax_transform(rank = "unique", trans = "clr") %>%
  ord_calc(
    method = "PCA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:10,
    colour="Facility", fill = "Facility",
    shape = "Year", alpha = 1,
    size = 3
  ) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=0.1))+
  theme(panel.grid = element_blank())+
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))
ggsave("PCA_Bacteria_loadings.png", plot =PC_16s_loadings, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Bacteria_loadings.svg", plot =PC_16s_loadings, device="svg", width=6, height=5, units="in",dpi=600)


PC_ITS_loadings<-phyloseq_ITS_clean %>%
  tax_transform(rank = "unique", trans = "clr") %>%
  ord_calc(
    method = "PCA"
  ) %>% 
  ord_plot(
    axes = c(1, 2),
    plot_taxa = 1:10,
    colour="Facility", fill = "Facility",
    shape = "Year", alpha = 1,
    size = 3
  ) + 
  scale_fill_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA),
        panel.border = element_rect(color='black',fill = NA,size=0.1))+
  theme(panel.grid = element_blank())+
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))
ggsave("PCA_Fungi_loadings.png", plot =PC_ITS_loadings, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Fungi_loadings.svg", plot =PC_ITS_loadings, device="svg", width=6, height=5, units="in",dpi=600)

##Note: will remove extra unlabeled loadings in Illustrator


# PERMANOVA #
#Calculate Aitchinson distance
dist_16s<-dist(asv.n0.clr_16s, method='euclidean')
dist_ITS<-dist(asv.n0.clr_ITS, method='euclidean')


#16s
permanova_16s<-pairwise.adonis2(dist_16s~Facility+Year+Facility:Year, data=metadata_16s_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_16s

#ITS
permanova_ITS<-pairwise.adonis2(dist_ITS~Facility+Year+Facility:Year, data=metadata_ITS_clean, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITS


##NOTE: Permanova showed significant difference by facility, season, and their interaction.

##Continue by (1) splitting ASV table by Facility effect within each facility and year 

#### Composition of microbiota - Barplots #### 
#ASV level plots

#Transform sample counts into compositions
asv_n0.acomp_16s<-as.data.frame(acomp(t(asv_16s_clean), total=1))
asv_n0.acomp_ITS<-as.data.frame(acomp(t(asv_ITS_clean), total=1))


#Make Phyloseq object for both years
phyloseq_16s_acomp<-phyloseq(otu_table(asv_n0.acomp_16s, taxa_are_rows = FALSE), tax_table(taxon_16s), sample_data(metadata_16s_clean))
phyloseq_ITS_acomp<-phyloseq(otu_table(asv_n0.acomp_ITS, taxa_are_rows = FALSE), tax_table(taxon_ITS), sample_data(metadata_ITS_clean))


#Make long format table from Phyloseq object combining all ASVs to the genus level
asv_16s_long <- phyloseq_16s_acomp %>%  
  transform_sample_counts(function(x) {x * 100} ) %>%
  tax_glom("Genus")%>%
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))

asv_ITS_long <- phyloseq_ITS_acomp %>%  
  transform_sample_counts(function(x) {x * 100} ) %>%
  tax_glom("Genus")%>% 
  psmelt() %>%  #Melts to long format
  arrange(desc(Abundance))


# #Make vector with ASV names above 1 or 10% relative abundance in at least one sample
# asv_16s_over5<-unique(c(asv_16s_long$OTU[which(asv_16s_long$Abundance >=8)]))
# asv_ITS_over5<-unique(c(asv_ITS_long$OTU[which(asv_ITS_long$Abundance >=8)])) 
# 
# #Filter table to obtain only ASVs with over 1 or 10% in at least one sample
# asv_16s_over5abund <- filter(asv_16s_long, OTU %in% asv_16s_over5)
# asv_ITS_over5abund <- filter(asv_ITS_long, OTU %in% asv_ITS_over5)

### Calculate mean relative abundance by Facility each ASV
asv_long_16s_mean<-asv_16s_long%>%
  group_by(Facility, Genus, Treatment, Time, Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_long_ITS_mean<-asv_ITS_long%>%
  group_by(Facility, Genus, Treatment, Time, Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

asv_long_16s_over5<-subset(asv_long_16s_mean, Mean>4)
asv_long_ITS_over5<-subset(asv_long_ITS_mean, Mean>8)


#Barplots
#Note: select distinct colors using: https://medialab.github.io/iwanthue/

barplot_16s<-ggplot(asv_long_16s_over5, aes(x=Time, y=Mean, fill=Genus))+
  geom_bar(stat='identity',color='black')+
  facet_grid(Year~Facility+Treatment, scales = "free", space = 'free')+
  ylab("Relative abundance (%)")+
  scale_x_discrete(limits=rev)+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10),
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 10)) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  scale_y_continuous(limits = c(0,100))+
  scale_fill_manual(values=c("#55bc4c","#e27a5b","#bd54c2","#8bb939","#7369de","#b3b13c","#8054a6","#d49f2b","#5e74c0","#db6e29",
                             "#3abbcc","#cc3c3a","#54bd7b","#d84589","#4a8532","#d288ca","#777522","#769eda",
                             "#53c2a3","#a64455","#6d733b","#36815b","#9eb56c","#9a4d7d","#e17e89","#d39f62","#995e2b"))

barplot_16s
ggsave("barplot_16s.png", plot=barplot_16s, device="png", width=8, height=10, units="in", dpi=600)
ggsave("barplot_16s.svg", plot=barplot_16s, device="svg", width=8, height=10, units="in", dpi=600)

barplot_ITS<-ggplot(asv_long_ITS_over5, aes(x=Time, y=Mean, fill=Genus))+
  geom_bar(stat='identity',color='black')+facet_grid(Year~Facility+Treatment, scales = "free", space = 'free')+
  ylab("Relative abundance (%)")+
  scale_x_discrete(limits=rev)+
  theme(axis.title.x = element_blank(), axis.text.x=element_text(color='black', size=10),
        axis.ticks=element_line(color='black'),
        axis.text.y = element_text(color='black', size=10)) +
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, nrow = 10)) +
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  #scale_color_viridis_d(begin = 0.2, end = 0.8, option='inferno')+
  #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
  scale_y_continuous(limits = c(0,100))+
  scale_fill_manual(values=c("#cc4437","#59b64c","#bc48ad","#a7b338","#705eca","#d19f45","#bd7dd8","#3a8250","#da3f83",
                             "#58c39b","#cb485e","#45aecf","#6d7127","#6d82c9","#db7030","#d681b2","#99af66","#9b486b",
                             "#9d6131","#de8777"))
  
barplot_ITS
ggsave("barplot_ITS.png", plot=barplot_ITS, device="png", width=8, height=10, units="in", dpi=600)
ggsave("barplot_ITS.svg", plot=barplot_ITS, device="svg", width=8, height=10, units="in", dpi=600)


#Calculate most abundant genera by facility, time and treatment.
#Calculate the mean relative abundance of genera by year and facility
#Make long format table from Phyloseq object
genus_16s<-asv_16s_long%>%
  group_by(Genus,Facility, Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))
write.csv(genus_16s, file="genus_16s.csv")

genus_ITS<-asv_ITS_long%>%
  group_by(Genus,Facility, Year)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))
write.csv(genus_ITS, file="genus_ITS.csv")

#Make vector with ASV names above 1 or 10% relative abundance in at least one sample
genus_16s_over5<-unique(c(genus_16s_long$OTU[which(genus_16s_long$Abundance >=5)]))
genus_ITS_over5<-unique(c(genus_ITS_long$OTU[which(genus_ITS_long$Abundance >=5)])) 

#Filter table to obtain only ASVs with over 1 or 10% in at least one sample
genus_16s_over5abund <- filter(genus_16s_long, OTU %in% genus_16s_over5)
genus_ITS_over5abund <- filter(genus_ITS_long, OTU %in% genus_ITS_over5)

### Calculate mean relative abundance by Facility each ASV
genus_16s_over5abund_mean<-genus_16s_over5abund%>%
  group_by(OTU, Facility, Genus, Treatment, Time)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))

genus_ITS_over5abund_mean<-genus_ITS_over5abund%>%
  group_by(OTU, Facility, Genus, Treatment, Time)%>%
  summarize(Mean=mean(Abundance))%>%
  arrange(desc(Mean))


# 
# heatmap_16s<-ggplot(genus_16s_over5abund_mean, aes(x=Treatment, y=Genus, fill=Mean))+
#   geom_tile(color='black')+facet_grid(.~Facility+Time)+
#   scale_size(limits=c(0,50),name = "Relative abundance (%)")+
#   theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10, angle=90), 
#         axis.ticks=element_line(color='black'),
#         axis.text.y = element_text(color='black', size=10)) + 
#   theme(legend.position = "bottom")+
#   theme(panel.background = element_rect(fill="grey99", color =NA),
#         plot.background = element_rect(fill="white", color =NA)) +
#   theme(strip.background= element_rect(fill="white", color = "black"), 
#         strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
#   ggtitle("Bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_viridis_c(begin = 1, end = 0, option='inferno')
# #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
# heatmap_16s
# ggsave("Heatmap_genus16s.svg", plot=heatmap_16s, device="svg", width=8, height=10, units="in", dpi=600)
# ggsave("Heatmap_genus16s.png", plot=heatmap_16s, device="png", width=8, height=10, units="in", dpi=600)
# 
# heatmap_ITS<-ggplot(genus_ITS_over5abund_mean, aes(x=Treatment, y=Genus, fill=Mean))+
#   geom_tile(color='black')+facet_grid(.~Facility+Time)+
#   scale_size(limits=c(0,50),name = "Relative abundance (%)")+
#   theme(axis.title = element_blank(), axis.text.x=element_text(color='black', size=10, angle=90), 
#         axis.ticks=element_line(color='black'),
#         axis.text.y = element_text(color='black', size=10)) + 
#   theme(legend.position = "bottom")+
#   theme(panel.background = element_rect(fill="grey99", color =NA),
#         plot.background = element_rect(fill="white", color =NA)) +
#   theme(strip.background= element_rect(fill="white", color = "black"), 
#         strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
#         panel.border = element_rect(color="black", fill=NA))+
#   theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
#   ggtitle("Fungi")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_viridis_c(begin = 1, end = 0, option='inferno')
# #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))
# heatmap_ITS
# ggsave("Heatmap_genusITS.svg", plot=heatmap_ITS, device="svg", width=8, height=10, units="in", dpi=600)
# ggsave("Heatmap_genusITS.png", plot=heatmap_ITS, device="png", width=8, height=10, units="in", dpi=600)


#Plots for selected groups of microorganisms
entero_16s<-subset(asv_16s_long, Family=="Enterobacteriaceae")
pseudo_16s<-subset(asv_16s_long, Genus=="Pseudomonas")
spore_16s<-subset(asv_16s_long, Genus=="Bacillus" | Genus =="Clostridium_sensu_stricto_1")

#Calculate summay statistics
entero_16s_mean<-entero_16s%>%
  group_by(Facility, Treatment, Time)%>%
  summarize(Mean=mean(Abundance), SD=sd(Abundance), SE=SD/sqrt(n()))%>%
  arrange(desc(Mean))

pseudo_16s_mean<-pseudo_16s%>%
  group_by(Facility, Treatment, Time)%>%
  summarize(Mean=mean(Abundance), SD=sd(Abundance), SE=SD/sqrt(n()))%>%
  arrange(desc(Mean))

spore_16s_mean<-spore_16s%>%
  group_by(Facility, Treatment, Time)%>%
  summarize(Mean=mean(Abundance), SD=sd(Abundance), SE=SD/sqrt(n()))%>%
  arrange(desc(Mean))

#plot
entero_plot<-ggplot(entero_16s_mean, aes(x=Time, y=Mean, fill=Treatment))+
  geom_bar(stat="identity", color="black")+geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2,position=position_dodge(.9)) +
  facet_grid(Facility~Treatment, scales = "free_y")+
  scale_x_discrete(limits=rev)+
  scale_fill_viridis_d(option='viridis')+
  theme(axis.text=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'))+
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Enterobacteriaceae")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  ylab("Mean relative abundance (%)")+ xlab("Sample collection time")


pseudo_plot<-ggplot(pseudo_16s_mean, aes(x=Time, y=Mean, fill=Treatment))+
  geom_bar(stat="identity", color="black")+geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2,position=position_dodge(.9)) +
  facet_grid(Facility~Treatment, scales = "free_y")+
  scale_x_discrete(limits=rev)+
  scale_fill_viridis_d(option='viridis')+
  theme(axis.text=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'))+
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Pseudomonas")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  ylab("Mean relative abundance (%)")+ xlab("Sample collection time")

spore_plot<-ggplot(spore_16s_mean, aes(x=Time, y=Mean, fill=Treatment))+
  geom_bar(stat="identity", color="black")+geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2,position=position_dodge(.9)) +
  facet_grid(Facility~Treatment, scales = "free_y")+
  scale_x_discrete(limits=rev)+
  scale_fill_viridis_d(option='viridis')+
  theme(axis.text=element_text(color='black', size=10), 
        axis.ticks=element_line(color='black'))+
  theme(legend.position = "bottom")+
  theme(panel.background = element_rect(fill="grey99", color =NA),
        plot.background = element_rect(fill="white", color =NA)) +
  theme(strip.background= element_rect(fill="white", color = "black"), 
        strip.text.y = element_text(size=10, angle=0, face = 'italic'), 
        panel.border = element_rect(color="black", fill=NA))+
  theme(legend.text=element_text(size=10,color='black'), legend.title= element_text(size=10)) +
  ggtitle("Spore forming bacteria")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
  ylab("Mean relative abundance (%)")+ xlab("Sample collection time")

ggsave("Enterobacteriaceae.png", plot=entero_plot, device="png", width=8, height=10, units="in", dpi=600)
ggsave("Pseudomonas.png", plot=pseudo_plot, device="png", width=8, height=10, units="in", dpi=600)
ggsave("Sporeformers.png", plot=spore_plot, device="png", width=8, height=10, units="in", dpi=600)





#### Analysis of microbiota composition by season, facility and treatment ####
#Filter samples by facility, Year, and Treatment
ps.16sF1T1Y1<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T1" & Year =="Y1")
ps.16sF1T2Y1<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T2" & Year =="Y1")
ps.16sF1T3Y1<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T3" & Year =="Y1")
ps.16sF1T4Y1<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T4" & Year =="Y1")

ps.16sF2T1Y1<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T1" & Year =="Y1")
ps.16sF2T2Y1<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T2" & Year =="Y1")
ps.16sF2T3Y1<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T3" & Year =="Y1")
ps.16sF2T4Y1<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T4" & Year =="Y1")

ps.16sF3T1Y1<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T1" & Year =="Y1")
ps.16sF3T2Y1<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T2" & Year =="Y1")
ps.16sF3T3Y1<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T3" & Year =="Y1")
ps.16sF3T4Y1<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T4" & Year =="Y1")

ps.16sF1T1Y2<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T1" & Year =="Y2")
ps.16sF1T2Y2<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T2" & Year =="Y2")
ps.16sF1T3Y2<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T3" & Year =="Y2")
ps.16sF1T4Y2<-subset_samples(phyloseq_16s_clean, Facility=="P1" & Treatment =="T4" & Year =="Y2")

ps.16sF2T1Y2<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T1" & Year =="Y2")
ps.16sF2T2Y2<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T2" & Year =="Y2")
ps.16sF2T3Y2<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T3" & Year =="Y2")
ps.16sF2T4Y2<-subset_samples(phyloseq_16s_clean, Facility=="P2" & Treatment =="T4" & Year =="Y2")

ps.16sF3T1Y2<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T1" & Year =="Y2")
ps.16sF3T2Y2<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T2" & Year =="Y2")
ps.16sF3T3Y2<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T3" & Year =="Y2")
ps.16sF3T4Y2<-subset_samples(phyloseq_16s_clean, Facility=="P3" & Treatment =="T4" & Year =="Y2")

ps.ITSF1T1Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T1" & Year =="Y1")
ps.ITSF1T2Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T2" & Year =="Y1")
ps.ITSF1T3Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T3" & Year =="Y1")
ps.ITSF1T4Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T4" & Year =="Y1")

ps.ITSF2T1Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T1" & Year =="Y1")
ps.ITSF2T2Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T2" & Year =="Y1")
ps.ITSF2T3Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T3" & Year =="Y1")
ps.ITSF2T4Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T4" & Year =="Y1")

ps.ITSF3T1Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T1" & Year =="Y1")
ps.ITSF3T2Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T2" & Year =="Y1")
ps.ITSF3T3Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T3" & Year =="Y1")
ps.ITSF3T4Y1<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T4" & Year =="Y1")

ps.ITSF1T1Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T1" & Year =="Y2")
ps.ITSF1T2Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T2" & Year =="Y2")
ps.ITSF1T3Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T3" & Year =="Y2")
ps.ITSF1T4Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P1" & Treatment =="T4" & Year =="Y2")

ps.ITSF2T1Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T1" & Year =="Y2")
ps.ITSF2T2Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T2" & Year =="Y2")
ps.ITSF2T3Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T3" & Year =="Y2")
ps.ITSF2T4Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P2" & Treatment =="T4" & Year =="Y2")

ps.ITSF3T1Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T1" & Year =="Y2")
ps.ITSF3T2Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T2" & Year =="Y2")
ps.ITSF3T3Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T3" & Year =="Y2")
ps.ITSF3T4Y2<-subset_samples(phyloseq_ITS_clean, Facility=="P3" & Treatment =="T4" & Year =="Y2")

#Get ASV table from phyloseq object
asv_16sF1T1Y1<-as.data.frame(t(otu_table(ps.16sF1T1Y1)))
asv_16sF1T2Y1<-as.data.frame(t(otu_table(ps.16sF1T2Y1)))
asv_16sF1T3Y1<-as.data.frame(t(otu_table(ps.16sF1T3Y1)))
asv_16sF1T4Y1<-as.data.frame(t(otu_table(ps.16sF1T4Y1)))

asv_16sF2T1Y1<-as.data.frame(t(otu_table(ps.16sF2T1Y1)))
asv_16sF2T2Y1<-as.data.frame(t(otu_table(ps.16sF2T2Y1)))
asv_16sF2T3Y1<-as.data.frame(t(otu_table(ps.16sF2T3Y1)))
asv_16sF2T4Y1<-as.data.frame(t(otu_table(ps.16sF2T4Y1)))

asv_16sF3T1Y1<-as.data.frame(t(otu_table(ps.16sF3T1Y1)))
asv_16sF3T2Y1<-as.data.frame(t(otu_table(ps.16sF3T2Y1)))
asv_16sF3T3Y1<-as.data.frame(t(otu_table(ps.16sF3T3Y1)))
asv_16sF3T4Y1<-as.data.frame(t(otu_table(ps.16sF3T4Y1)))

asv_16sF1T1Y2<-as.data.frame(t(otu_table(ps.16sF1T1Y2)))
asv_16sF1T2Y2<-as.data.frame(t(otu_table(ps.16sF1T2Y2)))
asv_16sF1T3Y2<-as.data.frame(t(otu_table(ps.16sF1T3Y2)))
asv_16sF1T4Y2<-as.data.frame(t(otu_table(ps.16sF1T4Y2)))

asv_16sF2T1Y2<-as.data.frame(t(otu_table(ps.16sF2T1Y2)))
asv_16sF2T2Y2<-as.data.frame(t(otu_table(ps.16sF2T2Y2)))
asv_16sF2T3Y2<-as.data.frame(t(otu_table(ps.16sF2T3Y2)))
asv_16sF2T4Y2<-as.data.frame(t(otu_table(ps.16sF2T4Y2)))

asv_16sF3T1Y2<-as.data.frame(t(otu_table(ps.16sF3T1Y2)))
asv_16sF3T2Y2<-as.data.frame(t(otu_table(ps.16sF3T2Y2)))
asv_16sF3T3Y2<-as.data.frame(t(otu_table(ps.16sF3T3Y2)))
asv_16sF3T4Y2<-as.data.frame(t(otu_table(ps.16sF3T4Y2)))

asv_ITSF1T1Y1<-as.data.frame(t(otu_table(ps.ITSF1T1Y1)))
asv_ITSF1T2Y1<-as.data.frame(t(otu_table(ps.ITSF1T2Y1)))
asv_ITSF1T3Y1<-as.data.frame(t(otu_table(ps.ITSF1T3Y1)))
asv_ITSF1T4Y1<-as.data.frame(t(otu_table(ps.ITSF1T4Y1)))

asv_ITSF2T1Y1<-as.data.frame(t(otu_table(ps.ITSF2T1Y1)))
asv_ITSF2T2Y1<-as.data.frame(t(otu_table(ps.ITSF2T2Y1)))
asv_ITSF2T3Y1<-as.data.frame(t(otu_table(ps.ITSF2T3Y1)))
asv_ITSF2T4Y1<-as.data.frame(t(otu_table(ps.ITSF2T4Y1)))

asv_ITSF3T1Y1<-as.data.frame(t(otu_table(ps.ITSF3T1Y1)))
asv_ITSF3T2Y1<-as.data.frame(t(otu_table(ps.ITSF3T2Y1)))
asv_ITSF3T3Y1<-as.data.frame(t(otu_table(ps.ITSF3T3Y1)))
asv_ITSF3T4Y1<-as.data.frame(t(otu_table(ps.ITSF3T4Y1)))

asv_ITSF1T1Y2<-as.data.frame(t(otu_table(ps.ITSF1T1Y2)))
asv_ITSF1T2Y2<-as.data.frame(t(otu_table(ps.ITSF1T2Y2)))
asv_ITSF1T3Y2<-as.data.frame(t(otu_table(ps.ITSF1T3Y2)))
asv_ITSF1T4Y2<-as.data.frame(t(otu_table(ps.ITSF1T4Y2)))

asv_ITSF2T1Y2<-as.data.frame(t(otu_table(ps.ITSF2T1Y2)))
asv_ITSF2T2Y2<-as.data.frame(t(otu_table(ps.ITSF2T2Y2)))
asv_ITSF2T3Y2<-as.data.frame(t(otu_table(ps.ITSF2T3Y2)))
asv_ITSF2T4Y2<-as.data.frame(t(otu_table(ps.ITSF2T4Y2)))

asv_ITSF3T1Y2<-as.data.frame(t(otu_table(ps.ITSF3T1Y2)))
asv_ITSF3T2Y2<-as.data.frame(t(otu_table(ps.ITSF3T2Y2)))
asv_ITSF3T3Y2<-as.data.frame(t(otu_table(ps.ITSF3T3Y2)))
asv_ITSF3T4Y2<-as.data.frame(t(otu_table(ps.ITSF3T4Y2)))

#Remove ASVs with zero counts in all samples
asv_16sF1T1Y1<-asv_16sF1T1Y1[ which(rowSums(asv_16sF1T1Y1)>0),]
asv_16sF1T1Y1<-t(asv_16sF1T1Y1)

asv_16sF1T2Y1<-asv_16sF1T2Y1[ which(rowSums(asv_16sF1T2Y1)>0),]
asv_16sF1T2Y1<-t(asv_16sF1T2Y1)

asv_16sF1T3Y1<-asv_16sF1T3Y1[ which(rowSums(asv_16sF1T3Y1)>0),]
asv_16sF1T3Y1<-t(asv_16sF1T3Y1)

asv_16sF1T4Y1<-asv_16sF1T4Y1[ which(rowSums(asv_16sF1T4Y1)>0),]
asv_16sF1T4Y1<-t(asv_16sF1T4Y1)

asv_16sF2T1Y1<-asv_16sF2T1Y1[ which(rowSums(asv_16sF2T1Y1)>0),]
asv_16sF2T1Y1<-t(asv_16sF2T1Y1)

asv_16sF2T2Y1<-asv_16sF2T2Y1[ which(rowSums(asv_16sF2T2Y1)>0),]
asv_16sF2T2Y1<-t(asv_16sF2T2Y1)

asv_16sF2T3Y1<-asv_16sF2T3Y1[ which(rowSums(asv_16sF2T3Y1)>0),]
asv_16sF2T3Y1<-t(asv_16sF2T3Y1)

asv_16sF2T4Y1<-asv_16sF2T4Y1[ which(rowSums(asv_16sF2T4Y1)>0),]
asv_16sF2T4Y1<-t(asv_16sF2T4Y1)

asv_16sF3T1Y1<-asv_16sF3T1Y1[ which(rowSums(asv_16sF3T1Y1)>0),]
asv_16sF3T1Y1<-t(asv_16sF3T1Y1)

asv_16sF3T2Y1<-asv_16sF3T2Y1[ which(rowSums(asv_16sF3T2Y1)>0),]
asv_16sF3T2Y1<-t(asv_16sF3T2Y1)

asv_16sF3T3Y1<-asv_16sF3T3Y1[ which(rowSums(asv_16sF3T3Y1)>0),]
asv_16sF3T3Y1<-t(asv_16sF3T3Y1)

asv_16sF3T4Y1<-asv_16sF3T4Y1[ which(rowSums(asv_16sF3T4Y1)>0),]
asv_16sF3T4Y1<-t(asv_16sF3T4Y1)

asv_16sF1T1Y2<-asv_16sF1T1Y2[ which(rowSums(asv_16sF1T1Y2)>0),]
asv_16sF1T1Y2<-t(asv_16sF1T1Y2)

asv_16sF1T2Y2<-asv_16sF1T2Y2[ which(rowSums(asv_16sF1T2Y2)>0),]
asv_16sF1T2Y2<-t(asv_16sF1T2Y2)

asv_16sF1T3Y2<-asv_16sF1T3Y2[ which(rowSums(asv_16sF1T3Y2)>0),]
asv_16sF1T3Y2<-t(asv_16sF1T3Y2)

asv_16sF1T4Y2<-asv_16sF1T4Y2[ which(rowSums(asv_16sF1T4Y2)>0),]
asv_16sF1T4Y2<-t(asv_16sF1T4Y2)

asv_16sF2T1Y2<-asv_16sF2T1Y2[ which(rowSums(asv_16sF2T1Y2)>0),]
asv_16sF2T1Y2<-t(asv_16sF2T1Y2)

asv_16sF2T2Y2<-asv_16sF2T2Y2[ which(rowSums(asv_16sF2T2Y2)>0),]
asv_16sF2T2Y2<-t(asv_16sF2T2Y2)

asv_16sF2T3Y2<-asv_16sF2T3Y2[ which(rowSums(asv_16sF2T3Y2)>0),]
asv_16sF2T3Y2<-t(asv_16sF2T3Y2)

asv_16sF2T4Y2<-asv_16sF2T4Y2[ which(rowSums(asv_16sF2T4Y2)>0),]
asv_16sF2T4Y2<-t(asv_16sF2T4Y2)

asv_16sF3T1Y2<-asv_16sF3T1Y2[ which(rowSums(asv_16sF3T1Y2)>0),]
asv_16sF3T1Y2<-t(asv_16sF3T1Y2)

asv_16sF3T2Y2<-asv_16sF3T2Y2[ which(rowSums(asv_16sF3T2Y2)>0),]
asv_16sF3T2Y2<-t(asv_16sF3T2Y2)

asv_16sF3T3Y2<-asv_16sF3T3Y2[ which(rowSums(asv_16sF3T3Y2)>0),]
asv_16sF3T3Y2<-t(asv_16sF3T3Y2)

asv_16sF3T4Y2<-asv_16sF3T4Y2[ which(rowSums(asv_16sF3T4Y2)>0),]
asv_16sF3T4Y2<-t(asv_16sF3T4Y2)

asv_ITSF1T1Y1<-asv_ITSF1T1Y1[ which(rowSums(asv_ITSF1T1Y1)>0),]
asv_ITSF1T1Y1<-t(asv_ITSF1T1Y1)

asv_ITSF1T2Y1<-asv_ITSF1T2Y1[ which(rowSums(asv_ITSF1T2Y1)>0),]
asv_ITSF1T2Y1<-t(asv_ITSF1T2Y1)

asv_ITSF1T3Y1<-asv_ITSF1T3Y1[ which(rowSums(asv_ITSF1T3Y1)>0),]
asv_ITSF1T3Y1<-t(asv_ITSF1T3Y1)

asv_ITSF1T4Y1<-asv_ITSF1T4Y1[ which(rowSums(asv_ITSF1T4Y1)>0),]
asv_ITSF1T4Y1<-t(asv_ITSF1T4Y1)

asv_ITSF2T1Y1<-asv_ITSF2T1Y1[ which(rowSums(asv_ITSF2T1Y1)>0),]
asv_ITSF2T1Y1<-t(asv_ITSF2T1Y1)

asv_ITSF2T2Y1<-asv_ITSF2T2Y1[ which(rowSums(asv_ITSF2T2Y1)>0),]
asv_ITSF2T2Y1<-t(asv_ITSF2T2Y1)

asv_ITSF2T3Y1<-asv_ITSF2T3Y1[ which(rowSums(asv_ITSF2T3Y1)>0),]
asv_ITSF2T3Y1<-t(asv_ITSF2T3Y1)

asv_ITSF2T4Y1<-asv_ITSF2T4Y1[ which(rowSums(asv_ITSF2T4Y1)>0),]
asv_ITSF2T4Y1<-t(asv_ITSF2T4Y1)

asv_ITSF3T1Y1<-asv_ITSF3T1Y1[ which(rowSums(asv_ITSF3T1Y1)>0),]
asv_ITSF3T1Y1<-t(asv_ITSF3T1Y1)

asv_ITSF3T2Y1<-asv_ITSF3T2Y1[ which(rowSums(asv_ITSF3T2Y1)>0),]
asv_ITSF3T2Y1<-t(asv_ITSF3T2Y1)

asv_ITSF3T3Y1<-asv_ITSF3T3Y1[ which(rowSums(asv_ITSF3T3Y1)>0),]
asv_ITSF3T3Y1<-t(asv_ITSF3T3Y1)

asv_ITSF3T4Y1<-asv_ITSF3T4Y1[ which(rowSums(asv_ITSF3T4Y1)>0),]
asv_ITSF3T4Y1<-t(asv_ITSF3T4Y1)

asv_ITSF1T1Y2<-asv_ITSF1T1Y2[ which(rowSums(asv_ITSF1T1Y2)>0),]
asv_ITSF1T1Y2<-t(asv_ITSF1T1Y2)

asv_ITSF1T2Y2<-asv_ITSF1T2Y2[ which(rowSums(asv_ITSF1T2Y2)>0),]
asv_ITSF1T2Y2<-t(asv_ITSF1T2Y2)

asv_ITSF1T3Y2<-asv_ITSF1T3Y2[ which(rowSums(asv_ITSF1T3Y2)>0),]
asv_ITSF1T3Y2<-t(asv_ITSF1T3Y2)

asv_ITSF1T4Y2<-asv_ITSF1T4Y2[ which(rowSums(asv_ITSF1T4Y2)>0),]
asv_ITSF1T4Y2<-t(asv_ITSF1T4Y2)

asv_ITSF2T1Y2<-asv_ITSF2T1Y2[ which(rowSums(asv_ITSF2T1Y2)>0),]
asv_ITSF2T1Y2<-t(asv_ITSF2T1Y2)

asv_ITSF2T2Y2<-asv_ITSF2T2Y2[ which(rowSums(asv_ITSF2T2Y2)>0),]
asv_ITSF2T2Y2<-t(asv_ITSF2T2Y2)

asv_ITSF2T3Y2<-asv_ITSF2T3Y2[ which(rowSums(asv_ITSF2T3Y2)>0),]
asv_ITSF2T3Y2<-t(asv_ITSF2T3Y2)

asv_ITSF2T4Y2<-asv_ITSF2T4Y2[ which(rowSums(asv_ITSF2T4Y2)>0),]
asv_ITSF2T4Y2<-t(asv_ITSF2T4Y2)

asv_ITSF3T1Y2<-asv_ITSF3T1Y2[ which(rowSums(asv_ITSF3T1Y2)>0),]
asv_ITSF3T1Y2<-t(asv_ITSF3T1Y2)

asv_ITSF3T2Y2<-asv_ITSF3T2Y2[ which(rowSums(asv_ITSF3T2Y2)>0),]
asv_ITSF3T2Y2<-t(asv_ITSF3T2Y2)

asv_ITSF3T3Y2<-asv_ITSF3T3Y2[ which(rowSums(asv_ITSF3T3Y2)>0),]
asv_ITSF3T3Y2<-t(asv_ITSF3T3Y2)

asv_ITSF3T4Y2<-asv_ITSF3T4Y2[ which(rowSums(asv_ITSF3T4Y2)>0),]
asv_ITSF3T4Y2<-t(asv_ITSF3T4Y2)

#Get metadata 
metadata.16sF1T1Y1<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T1" & Year =="Y1")
metadata.16sF1T2Y1<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T2" & Year =="Y1")
metadata.16sF1T3Y1<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T3" & Year =="Y1")
metadata.16sF1T4Y1<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T4" & Year =="Y1")

metadata.16sF2T1Y1<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T1" & Year =="Y1")
metadata.16sF2T2Y1<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T2" & Year =="Y1")
metadata.16sF2T3Y1<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T3" & Year =="Y1")
metadata.16sF2T4Y1<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T4" & Year =="Y1")

metadata.16sF3T1Y1<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T1" & Year =="Y1")
metadata.16sF3T2Y1<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T2" & Year =="Y1")
metadata.16sF3T3Y1<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T3" & Year =="Y1")
metadata.16sF3T4Y1<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T4" & Year =="Y1")

metadata.16sF1T1Y2<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T1" & Year =="Y2")
metadata.16sF1T2Y2<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T2" & Year =="Y2")
metadata.16sF1T3Y2<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T3" & Year =="Y2")
metadata.16sF1T4Y2<-subset(metadata_16s_clean, Facility == "P1" & Treatment =="T4" & Year =="Y2")

metadata.16sF2T1Y2<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T1" & Year =="Y2")
metadata.16sF2T2Y2<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T2" & Year =="Y2")
metadata.16sF2T3Y2<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T3" & Year =="Y2")
metadata.16sF2T4Y2<-subset(metadata_16s_clean, Facility == "P2" & Treatment =="T4" & Year =="Y2")

metadata.16sF3T1Y2<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T1" & Year =="Y2")
metadata.16sF3T2Y2<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T2" & Year =="Y2")
metadata.16sF3T3Y2<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T3" & Year =="Y2")
metadata.16sF3T4Y2<-subset(metadata_16s_clean, Facility == "P3" & Treatment =="T4" & Year =="Y2")

metadata.ITSF1T1Y1<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T1" & Year =="Y1")
metadata.ITSF1T2Y1<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T2" & Year =="Y1")
metadata.ITSF1T3Y1<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T3" & Year =="Y1")
metadata.ITSF1T4Y1<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T4" & Year =="Y1")

metadata.ITSF2T1Y1<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T1" & Year =="Y1")
metadata.ITSF2T2Y1<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T2" & Year =="Y1")
metadata.ITSF2T3Y1<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T3" & Year =="Y1")
metadata.ITSF2T4Y1<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T4" & Year =="Y1")

metadata.ITSF3T1Y1<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T1" & Year =="Y1")
metadata.ITSF3T2Y1<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T2" & Year =="Y1")
metadata.ITSF3T3Y1<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T3" & Year =="Y1")
metadata.ITSF3T4Y1<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T4" & Year =="Y1")

metadata.ITSF1T1Y2<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T1" & Year =="Y2")
metadata.ITSF1T2Y2<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T2" & Year =="Y2")
metadata.ITSF1T3Y2<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T3" & Year =="Y2")
metadata.ITSF1T4Y2<-subset(metadata_ITS_clean, Facility == "P1" & Treatment =="T4" & Year =="Y2")

metadata.ITSF2T1Y2<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T1" & Year =="Y2")
metadata.ITSF2T2Y2<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T2" & Year =="Y2")
metadata.ITSF2T3Y2<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T3" & Year =="Y2")
metadata.ITSF2T4Y2<-subset(metadata_ITS_clean, Facility == "P2" & Treatment =="T4" & Year =="Y2")

metadata.ITSF3T1Y2<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T1" & Year =="Y2")
metadata.ITSF3T2Y2<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T2" & Year =="Y2")
metadata.ITSF3T3Y2<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T3" & Year =="Y2")
metadata.ITSF3T4Y2<-subset(metadata_ITS_clean, Facility == "P3" & Treatment =="T4" & Year =="Y2")


#Step 1: Convert ASV table to appropriate format
#Following step requires samples on rows and ASVs in columns
head(asv_16sF1T1Y1) 
head(asv_16sF1T2Y1) 
head(asv_16sF1T3Y1) 
head(asv_16sF1T4Y1) 

head(asv_16sF2T1Y1) 
head(asv_16sF2T2Y1) 
head(asv_16sF2T3Y1) 
head(asv_16sF2T4Y1) 

head(asv_16sF3T1Y1) 
head(asv_16sF3T2Y1) 
head(asv_16sF3T3Y1) 
head(asv_16sF3T4Y1) 

head(asv_16sF1T1Y2) 
head(asv_16sF1T2Y2) 
head(asv_16sF1T3Y2) 
head(asv_16sF1T4Y2) 

head(asv_16sF2T1Y2) 
head(asv_16sF2T2Y2) 
head(asv_16sF2T3Y2) 
head(asv_16sF2T4Y2) 

head(asv_16sF3T1Y2) 
head(asv_16sF3T2Y2) 
head(asv_16sF3T3Y2) 
head(asv_16sF3T4Y2) 

head(asv_ITSF1T1Y1) 
head(asv_ITSF1T2Y1) 
head(asv_ITSF1T3Y1) 
head(asv_ITSF1T4Y1) 

head(asv_ITSF2T1Y1) 
head(asv_ITSF2T2Y1) 
head(asv_ITSF2T3Y1) 
head(asv_ITSF2T4Y1) 

head(asv_ITSF3T1Y1) 
head(asv_ITSF3T2Y1) 
head(asv_ITSF3T3Y1) 
head(asv_ITSF3T4Y1) 

head(asv_ITSF1T1Y2) 
head(asv_ITSF1T2Y2) 
head(asv_ITSF1T3Y2) 
head(asv_ITSF1T4Y2) 

head(asv_ITSF2T1Y2) 
head(asv_ITSF2T2Y2) 
head(asv_ITSF2T3Y2) 
head(asv_ITSF2T4Y2) 

head(asv_ITSF3T1Y2) 
head(asv_ITSF3T2Y2) 
head(asv_ITSF3T3Y2) 
head(asv_ITSF3T4Y2)

#Step 2: Replace zero values before clr transformation. 
#Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
asv.n0_16sF1T1Y1<-t(cmultRepl(asv_16sF1T1Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF1T2Y1<-t(cmultRepl(asv_16sF1T2Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF1T3Y1<-t(cmultRepl(asv_16sF1T3Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF1T4Y1<-t(cmultRepl(asv_16sF1T4Y1, label=0, method="CZM", output="p-counts")) 

asv.n0_16sF2T1Y1<-t(cmultRepl(asv_16sF2T1Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF2T2Y1<-t(cmultRepl(asv_16sF2T2Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF2T3Y1<-t(cmultRepl(asv_16sF2T3Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF2T4Y1<-t(cmultRepl(asv_16sF2T4Y1, label=0, method="CZM", output="p-counts")) 

asv.n0_16sF3T1Y1<-t(cmultRepl(asv_16sF3T1Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF3T2Y1<-t(cmultRepl(asv_16sF3T2Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF3T3Y1<-t(cmultRepl(asv_16sF3T3Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF3T4Y1<-t(cmultRepl(asv_16sF3T4Y1, label=0, method="CZM", output="p-counts")) 

asv.n0_16sF1T1Y2<-t(cmultRepl(asv_16sF1T1Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF1T2Y2<-t(cmultRepl(asv_16sF1T2Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF1T3Y2<-t(cmultRepl(asv_16sF1T3Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF1T4Y2<-t(cmultRepl(asv_16sF1T4Y2, label=0, method="CZM", output="p-counts")) 

asv.n0_16sF2T1Y2<-t(cmultRepl(asv_16sF2T1Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF2T2Y2<-t(cmultRepl(asv_16sF2T2Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF2T3Y2<-t(cmultRepl(asv_16sF2T3Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF2T4Y2<-t(cmultRepl(asv_16sF2T4Y2, label=0, method="CZM", output="p-counts")) 

asv.n0_16sF3T1Y2<-t(cmultRepl(asv_16sF3T1Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF3T2Y2<-t(cmultRepl(asv_16sF3T2Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF3T3Y2<-t(cmultRepl(asv_16sF3T3Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_16sF3T4Y2<-t(cmultRepl(asv_16sF3T4Y2, label=0, method="CZM", output="p-counts")) 

asv.n0_ITSF1T1Y1<-t(cmultRepl(asv_ITSF1T1Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF1T2Y1<-t(cmultRepl(asv_ITSF1T2Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF1T3Y1<-t(cmultRepl(asv_ITSF1T3Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF1T4Y1<-t(cmultRepl(asv_ITSF1T4Y1, label=0, method="CZM", output="p-counts")) 

asv.n0_ITSF2T1Y1<-t(cmultRepl(asv_ITSF2T1Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF2T2Y1<-t(cmultRepl(asv_ITSF2T2Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF2T3Y1<-t(cmultRepl(asv_ITSF2T3Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF2T4Y1<-t(cmultRepl(asv_ITSF2T4Y1, label=0, method="CZM", output="p-counts")) 

asv.n0_ITSF3T1Y1<-t(cmultRepl(asv_ITSF3T1Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF3T2Y1<-t(cmultRepl(asv_ITSF3T2Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF3T3Y1<-t(cmultRepl(asv_ITSF3T3Y1, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF3T4Y1<-t(cmultRepl(asv_ITSF3T4Y1, label=0, method="CZM", output="p-counts")) 

asv.n0_ITSF1T1Y2<-t(cmultRepl(asv_ITSF1T1Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF1T2Y2<-t(cmultRepl(asv_ITSF1T2Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF1T3Y2<-t(cmultRepl(asv_ITSF1T3Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF1T4Y2<-t(cmultRepl(asv_ITSF1T4Y2, label=0, method="CZM", output="p-counts")) 

asv.n0_ITSF2T1Y2<-t(cmultRepl(asv_ITSF2T1Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF2T2Y2<-t(cmultRepl(asv_ITSF2T2Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF2T3Y2<-t(cmultRepl(asv_ITSF2T3Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF2T4Y2<-t(cmultRepl(asv_ITSF2T4Y2, label=0, method="CZM", output="p-counts")) 

asv.n0_ITSF3T1Y2<-t(cmultRepl(asv_ITSF3T1Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF3T2Y2<-t(cmultRepl(asv_ITSF3T2Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF3T3Y2<-t(cmultRepl(asv_ITSF3T3Y2, label=0, method="CZM", output="p-counts")) 
asv.n0_ITSF3T4Y2<-t(cmultRepl(asv_ITSF3T4Y2, label=0, method="CZM", output="p-counts")) 

#Note: Check the output to make sure there are no negative numbers. If samples or ASV are sparse, the CZM method will 
#add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
#into positives
asv.n0_16sF1T1Y1<-ifelse(asv.n0_16sF1T1Y1 < 0, asv.n0_16sF1T1Y1*(-1), asv.n0_16sF1T1Y1)
asv.n0_16sF1T2Y1<-ifelse(asv.n0_16sF1T2Y1 < 0, asv.n0_16sF1T2Y1*(-1), asv.n0_16sF1T2Y1)
asv.n0_16sF1T3Y1<-ifelse(asv.n0_16sF1T3Y1 < 0, asv.n0_16sF1T3Y1*(-1), asv.n0_16sF1T3Y1)
asv.n0_16sF1T4Y1<-ifelse(asv.n0_16sF1T4Y1 < 0, asv.n0_16sF1T4Y1*(-1), asv.n0_16sF1T4Y1)

asv.n0_16sF2T1Y1<-ifelse(asv.n0_16sF2T1Y1 < 0, asv.n0_16sF2T1Y1*(-1), asv.n0_16sF2T1Y1)
asv.n0_16sF2T2Y1<-ifelse(asv.n0_16sF2T2Y1 < 0, asv.n0_16sF2T2Y1*(-1), asv.n0_16sF2T2Y1)
asv.n0_16sF2T3Y1<-ifelse(asv.n0_16sF2T3Y1 < 0, asv.n0_16sF2T3Y1*(-1), asv.n0_16sF2T3Y1)
asv.n0_16sF2T4Y1<-ifelse(asv.n0_16sF2T4Y1 < 0, asv.n0_16sF2T4Y1*(-1), asv.n0_16sF2T4Y1)

asv.n0_16sF3T1Y1<-ifelse(asv.n0_16sF3T1Y1 < 0, asv.n0_16sF3T1Y1*(-1), asv.n0_16sF3T1Y1)
asv.n0_16sF3T2Y1<-ifelse(asv.n0_16sF3T2Y1 < 0, asv.n0_16sF3T2Y1*(-1), asv.n0_16sF3T2Y1)
asv.n0_16sF3T3Y1<-ifelse(asv.n0_16sF3T3Y1 < 0, asv.n0_16sF3T3Y1*(-1), asv.n0_16sF3T3Y1)
asv.n0_16sF3T4Y1<-ifelse(asv.n0_16sF3T4Y1 < 0, asv.n0_16sF3T4Y1*(-1), asv.n0_16sF3T4Y1)

asv.n0_16sF1T1Y2<-ifelse(asv.n0_16sF1T1Y2 < 0, asv.n0_16sF1T1Y2*(-1), asv.n0_16sF1T1Y2)
asv.n0_16sF1T2Y2<-ifelse(asv.n0_16sF1T2Y2 < 0, asv.n0_16sF1T2Y2*(-1), asv.n0_16sF1T2Y2)
asv.n0_16sF1T3Y2<-ifelse(asv.n0_16sF1T3Y2 < 0, asv.n0_16sF1T3Y2*(-1), asv.n0_16sF1T3Y2)
asv.n0_16sF1T4Y2<-ifelse(asv.n0_16sF1T4Y2 < 0, asv.n0_16sF1T4Y2*(-1), asv.n0_16sF1T4Y2)

asv.n0_16sF2T1Y2<-ifelse(asv.n0_16sF2T1Y2 < 0, asv.n0_16sF2T1Y2*(-1), asv.n0_16sF2T1Y2)
asv.n0_16sF2T2Y2<-ifelse(asv.n0_16sF2T2Y2 < 0, asv.n0_16sF2T2Y2*(-1), asv.n0_16sF2T2Y2)
asv.n0_16sF2T3Y2<-ifelse(asv.n0_16sF2T3Y2 < 0, asv.n0_16sF2T3Y2*(-1), asv.n0_16sF2T3Y2)
asv.n0_16sF2T4Y2<-ifelse(asv.n0_16sF2T4Y2 < 0, asv.n0_16sF2T4Y2*(-1), asv.n0_16sF2T4Y2)

asv.n0_16sF3T1Y2<-ifelse(asv.n0_16sF3T1Y2 < 0, asv.n0_16sF3T1Y2*(-1), asv.n0_16sF3T1Y2)
asv.n0_16sF3T2Y2<-ifelse(asv.n0_16sF3T2Y2 < 0, asv.n0_16sF3T2Y2*(-1), asv.n0_16sF3T2Y2)
asv.n0_16sF3T3Y2<-ifelse(asv.n0_16sF3T3Y2 < 0, asv.n0_16sF3T3Y2*(-1), asv.n0_16sF3T3Y2)
asv.n0_16sF3T4Y2<-ifelse(asv.n0_16sF3T4Y2 < 0, asv.n0_16sF3T4Y2*(-1), asv.n0_16sF3T4Y2)

asv.n0_ITSF1T1Y1<-ifelse(asv.n0_ITSF1T1Y1 < 0, asv.n0_ITSF1T1Y1*(-1), asv.n0_ITSF1T1Y1)
asv.n0_ITSF1T2Y1<-ifelse(asv.n0_ITSF1T2Y1 < 0, asv.n0_ITSF1T2Y1*(-1), asv.n0_ITSF1T2Y1)
asv.n0_ITSF1T3Y1<-ifelse(asv.n0_ITSF1T3Y1 < 0, asv.n0_ITSF1T3Y1*(-1), asv.n0_ITSF1T3Y1)
asv.n0_ITSF1T4Y1<-ifelse(asv.n0_ITSF1T4Y1 < 0, asv.n0_ITSF1T4Y1*(-1), asv.n0_ITSF1T4Y1)

asv.n0_ITSF2T1Y1<-ifelse(asv.n0_ITSF2T1Y1 < 0, asv.n0_ITSF2T1Y1*(-1), asv.n0_ITSF2T1Y1)
asv.n0_ITSF2T2Y1<-ifelse(asv.n0_ITSF2T2Y1 < 0, asv.n0_ITSF2T2Y1*(-1), asv.n0_ITSF2T2Y1)
asv.n0_ITSF2T3Y1<-ifelse(asv.n0_ITSF2T3Y1 < 0, asv.n0_ITSF2T3Y1*(-1), asv.n0_ITSF2T3Y1)
asv.n0_ITSF2T4Y1<-ifelse(asv.n0_ITSF2T4Y1 < 0, asv.n0_ITSF2T4Y1*(-1), asv.n0_ITSF2T4Y1)

asv.n0_ITSF3T1Y1<-ifelse(asv.n0_ITSF3T1Y1 < 0, asv.n0_ITSF3T1Y1*(-1), asv.n0_ITSF3T1Y1)
asv.n0_ITSF3T2Y1<-ifelse(asv.n0_ITSF3T2Y1 < 0, asv.n0_ITSF3T2Y1*(-1), asv.n0_ITSF3T2Y1)
asv.n0_ITSF3T3Y1<-ifelse(asv.n0_ITSF3T3Y1 < 0, asv.n0_ITSF3T3Y1*(-1), asv.n0_ITSF3T3Y1)
asv.n0_ITSF3T4Y1<-ifelse(asv.n0_ITSF3T4Y1 < 0, asv.n0_ITSF3T4Y1*(-1), asv.n0_ITSF3T4Y1)

asv.n0_ITSF1T1Y2<-ifelse(asv.n0_ITSF1T1Y2 < 0, asv.n0_ITSF1T1Y2*(-1), asv.n0_ITSF1T1Y2)
asv.n0_ITSF1T2Y2<-ifelse(asv.n0_ITSF1T2Y2 < 0, asv.n0_ITSF1T2Y2*(-1), asv.n0_ITSF1T2Y2)
asv.n0_ITSF1T3Y2<-ifelse(asv.n0_ITSF1T3Y2 < 0, asv.n0_ITSF1T3Y2*(-1), asv.n0_ITSF1T3Y2)
asv.n0_ITSF1T4Y2<-ifelse(asv.n0_ITSF1T4Y2 < 0, asv.n0_ITSF1T4Y2*(-1), asv.n0_ITSF1T4Y2)

asv.n0_ITSF2T1Y2<-ifelse(asv.n0_ITSF2T1Y2 < 0, asv.n0_ITSF2T1Y2*(-1), asv.n0_ITSF2T1Y2)
asv.n0_ITSF2T2Y2<-ifelse(asv.n0_ITSF2T2Y2 < 0, asv.n0_ITSF2T2Y2*(-1), asv.n0_ITSF2T2Y2)
asv.n0_ITSF2T3Y2<-ifelse(asv.n0_ITSF2T3Y2 < 0, asv.n0_ITSF2T3Y2*(-1), asv.n0_ITSF2T3Y2)
asv.n0_ITSF2T4Y2<-ifelse(asv.n0_ITSF2T4Y2 < 0, asv.n0_ITSF2T4Y2*(-1), asv.n0_ITSF2T4Y2)

asv.n0_ITSF3T1Y2<-ifelse(asv.n0_ITSF3T1Y2 < 0, asv.n0_ITSF3T1Y2*(-1), asv.n0_ITSF3T1Y2)
asv.n0_ITSF3T2Y2<-ifelse(asv.n0_ITSF3T2Y2 < 0, asv.n0_ITSF3T2Y2*(-1), asv.n0_ITSF3T2Y2)
asv.n0_ITSF3T3Y2<-ifelse(asv.n0_ITSF3T3Y2 < 0, asv.n0_ITSF3T3Y2*(-1), asv.n0_ITSF3T3Y2)
asv.n0_ITSF3T4Y2<-ifelse(asv.n0_ITSF3T4Y2 < 0, asv.n0_ITSF3T4Y2*(-1), asv.n0_ITSF3T4Y2)

#output table needs to have samples in columns and ASVs in rows
head(asv.n0_16sF1T1Y1) 
head(asv.n0_16sF1T2Y1) 
head(asv.n0_16sF1T3Y1) 
head(asv.n0_16sF1T4Y1)

head(asv.n0_16sF2T1Y1) 
head(asv.n0_16sF2T2Y1) 
head(asv.n0_16sF2T3Y1) 
head(asv.n0_16sF2T4Y1)

head(asv.n0_16sF3T1Y1) 
head(asv.n0_16sF3T2Y1) 
head(asv.n0_16sF3T3Y1) 
head(asv.n0_16sF3T4Y1)

head(asv.n0_16sF1T1Y2) 
head(asv.n0_16sF1T2Y2) 
head(asv.n0_16sF1T3Y2) 
head(asv.n0_16sF1T4Y2)

head(asv.n0_16sF2T1Y2) 
head(asv.n0_16sF2T2Y2) 
head(asv.n0_16sF2T3Y2) 
head(asv.n0_16sF2T4Y2)

head(asv.n0_16sF3T1Y2) 
head(asv.n0_16sF3T2Y2) 
head(asv.n0_16sF3T3Y2) 
head(asv.n0_16sF3T4Y2)

head(asv.n0_ITSF1T1Y1) 
head(asv.n0_ITSF1T2Y1) 
head(asv.n0_ITSF1T3Y1) 
head(asv.n0_ITSF1T4Y1)

head(asv.n0_ITSF2T1Y1) 
head(asv.n0_ITSF2T2Y1) 
head(asv.n0_ITSF2T3Y1) 
head(asv.n0_ITSF2T4Y1)

head(asv.n0_ITSF3T1Y1) 
head(asv.n0_ITSF3T2Y1) 
head(asv.n0_ITSF3T3Y1) 
head(asv.n0_ITSF3T4Y1)

head(asv.n0_ITSF1T1Y2) 
head(asv.n0_ITSF1T2Y2) 
head(asv.n0_ITSF1T3Y2) 
head(asv.n0_ITSF1T4Y2)

head(asv.n0_ITSF2T1Y2) 
head(asv.n0_ITSF2T2Y2) 
head(asv.n0_ITSF2T3Y2) 
head(asv.n0_ITSF2T4Y2)

head(asv.n0_ITSF3T1Y2) 
head(asv.n0_ITSF3T2Y2) 
head(asv.n0_ITSF3T3Y2) 
head(asv.n0_ITSF3T4Y2)

#Step 3: Convert data to proportions
asv.n0_16sF1T1Y1_prop<-apply(asv.n0_16sF1T1Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF1T2Y1_prop<-apply(asv.n0_16sF1T2Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF1T3Y1_prop<-apply(asv.n0_16sF1T3Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF1T4Y1_prop<-apply(asv.n0_16sF1T4Y1, 2, function(x) {x/sum(x)})

asv.n0_16sF2T1Y1_prop<-apply(asv.n0_16sF2T1Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF2T2Y1_prop<-apply(asv.n0_16sF2T2Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF2T3Y1_prop<-apply(asv.n0_16sF2T3Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF2T4Y1_prop<-apply(asv.n0_16sF2T4Y1, 2, function(x) {x/sum(x)})

asv.n0_16sF3T1Y1_prop<-apply(asv.n0_16sF3T1Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF3T2Y1_prop<-apply(asv.n0_16sF3T2Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF3T3Y1_prop<-apply(asv.n0_16sF3T3Y1, 2, function(x) {x/sum(x)})
asv.n0_16sF3T4Y1_prop<-apply(asv.n0_16sF3T4Y1, 2, function(x) {x/sum(x)})

asv.n0_16sF1T1Y2_prop<-apply(asv.n0_16sF1T1Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF1T2Y2_prop<-apply(asv.n0_16sF1T2Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF1T3Y2_prop<-apply(asv.n0_16sF1T3Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF1T4Y2_prop<-apply(asv.n0_16sF1T4Y2, 2, function(x) {x/sum(x)})

asv.n0_16sF2T1Y2_prop<-apply(asv.n0_16sF2T1Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF2T2Y2_prop<-apply(asv.n0_16sF2T2Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF2T3Y2_prop<-apply(asv.n0_16sF2T3Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF2T4Y2_prop<-apply(asv.n0_16sF2T4Y2, 2, function(x) {x/sum(x)})

asv.n0_16sF3T1Y2_prop<-apply(asv.n0_16sF3T1Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF3T2Y2_prop<-apply(asv.n0_16sF3T2Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF3T3Y2_prop<-apply(asv.n0_16sF3T3Y2, 2, function(x) {x/sum(x)})
asv.n0_16sF3T4Y2_prop<-apply(asv.n0_16sF3T4Y2, 2, function(x) {x/sum(x)})

asv.n0_ITSF1T1Y1_prop<-apply(asv.n0_ITSF1T1Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF1T2Y1_prop<-apply(asv.n0_ITSF1T2Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF1T3Y1_prop<-apply(asv.n0_ITSF1T3Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF1T4Y1_prop<-apply(asv.n0_ITSF1T4Y1, 2, function(x) {x/sum(x)})

asv.n0_ITSF2T1Y1_prop<-apply(asv.n0_ITSF2T1Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF2T2Y1_prop<-apply(asv.n0_ITSF2T2Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF2T3Y1_prop<-apply(asv.n0_ITSF2T3Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF2T4Y1_prop<-apply(asv.n0_ITSF2T4Y1, 2, function(x) {x/sum(x)})

asv.n0_ITSF3T1Y1_prop<-apply(asv.n0_ITSF3T1Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF3T2Y1_prop<-apply(asv.n0_ITSF3T2Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF3T3Y1_prop<-apply(asv.n0_ITSF3T3Y1, 2, function(x) {x/sum(x)})
asv.n0_ITSF3T4Y1_prop<-apply(asv.n0_ITSF3T4Y1, 2, function(x) {x/sum(x)})

asv.n0_ITSF1T1Y2_prop<-apply(asv.n0_ITSF1T1Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF1T2Y2_prop<-apply(asv.n0_ITSF1T2Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF1T3Y2_prop<-apply(asv.n0_ITSF1T3Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF1T4Y2_prop<-apply(asv.n0_ITSF1T4Y2, 2, function(x) {x/sum(x)})

asv.n0_ITSF2T1Y2_prop<-apply(asv.n0_ITSF2T1Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF2T2Y2_prop<-apply(asv.n0_ITSF2T2Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF2T3Y2_prop<-apply(asv.n0_ITSF2T3Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF2T4Y2_prop<-apply(asv.n0_ITSF2T4Y2, 2, function(x) {x/sum(x)})

asv.n0_ITSF3T1Y2_prop<-apply(asv.n0_ITSF3T1Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF3T2Y2_prop<-apply(asv.n0_ITSF3T2Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF3T3Y2_prop<-apply(asv.n0_ITSF3T3Y2, 2, function(x) {x/sum(x)})
asv.n0_ITSF3T4Y2_prop<-apply(asv.n0_ITSF3T4Y2, 2, function(x) {x/sum(x)})


#Step 4: Perform abundance and sample filtering and deal sparsity
#Filter the data to remove all taxa with less than 0.00001% abundance in any sample
asv.n0_16sF1T1Y1_prop_f<-asv.n0_16sF1T1Y1[apply(asv.n0_16sF1T1Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF1T2Y1_prop_f<-asv.n0_16sF1T2Y1[apply(asv.n0_16sF1T2Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF1T3Y1_prop_f<-asv.n0_16sF1T3Y1[apply(asv.n0_16sF1T3Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF1T4Y1_prop_f<-asv.n0_16sF1T4Y1[apply(asv.n0_16sF1T4Y1_prop, 1, min) > 0.000001, ]

asv.n0_16sF2T1Y1_prop_f<-asv.n0_16sF2T1Y1[apply(asv.n0_16sF2T1Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF2T2Y1_prop_f<-asv.n0_16sF2T2Y1[apply(asv.n0_16sF2T2Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF2T3Y1_prop_f<-asv.n0_16sF2T3Y1[apply(asv.n0_16sF2T3Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF2T4Y1_prop_f<-asv.n0_16sF2T4Y1[apply(asv.n0_16sF2T4Y1_prop, 1, min) > 0.000001, ]

asv.n0_16sF3T1Y1_prop_f<-asv.n0_16sF3T1Y1[apply(asv.n0_16sF3T1Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF3T2Y1_prop_f<-asv.n0_16sF3T2Y1[apply(asv.n0_16sF3T2Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF3T3Y1_prop_f<-asv.n0_16sF3T3Y1[apply(asv.n0_16sF3T3Y1_prop, 1, min) > 0.000001, ]
asv.n0_16sF3T4Y1_prop_f<-asv.n0_16sF3T4Y1[apply(asv.n0_16sF3T4Y1_prop, 1, min) > 0.000001, ]

asv.n0_16sF1T1Y2_prop_f<-asv.n0_16sF1T1Y2[apply(asv.n0_16sF1T1Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF1T2Y2_prop_f<-asv.n0_16sF1T2Y2[apply(asv.n0_16sF1T2Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF1T3Y2_prop_f<-asv.n0_16sF1T3Y2[apply(asv.n0_16sF1T3Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF1T4Y2_prop_f<-asv.n0_16sF1T4Y2[apply(asv.n0_16sF1T4Y2_prop, 1, min) > 0.000001, ]

asv.n0_16sF2T1Y2_prop_f<-asv.n0_16sF2T1Y2[apply(asv.n0_16sF2T1Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF2T2Y2_prop_f<-asv.n0_16sF2T2Y2[apply(asv.n0_16sF2T2Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF2T3Y2_prop_f<-asv.n0_16sF2T3Y2[apply(asv.n0_16sF2T3Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF2T4Y2_prop_f<-asv.n0_16sF2T4Y2[apply(asv.n0_16sF2T4Y2_prop, 1, min) > 0.000001, ]

asv.n0_16sF3T1Y2_prop_f<-asv.n0_16sF3T1Y2[apply(asv.n0_16sF3T1Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF3T2Y2_prop_f<-asv.n0_16sF3T2Y2[apply(asv.n0_16sF3T2Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF3T3Y2_prop_f<-asv.n0_16sF3T3Y2[apply(asv.n0_16sF3T3Y2_prop, 1, min) > 0.000001, ]
asv.n0_16sF3T4Y2_prop_f<-asv.n0_16sF3T4Y2[apply(asv.n0_16sF3T4Y2_prop, 1, min) > 0.000001, ]

asv.n0_ITSF1T1Y1_prop_f<-asv.n0_ITSF1T1Y1[apply(asv.n0_ITSF1T1Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF1T2Y1_prop_f<-asv.n0_ITSF1T2Y1[apply(asv.n0_ITSF1T2Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF1T3Y1_prop_f<-asv.n0_ITSF1T3Y1[apply(asv.n0_ITSF1T3Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF1T4Y1_prop_f<-asv.n0_ITSF1T4Y1[apply(asv.n0_ITSF1T4Y1_prop, 1, min) > 0.000001, ]

asv.n0_ITSF2T1Y1_prop_f<-asv.n0_ITSF2T1Y1[apply(asv.n0_ITSF2T1Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF2T2Y1_prop_f<-asv.n0_ITSF2T2Y1[apply(asv.n0_ITSF2T2Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF2T3Y1_prop_f<-asv.n0_ITSF2T3Y1[apply(asv.n0_ITSF2T3Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF2T4Y1_prop_f<-asv.n0_ITSF2T4Y1[apply(asv.n0_ITSF2T4Y1_prop, 1, min) > 0.000001, ]

asv.n0_ITSF3T1Y1_prop_f<-asv.n0_ITSF3T1Y1[apply(asv.n0_ITSF3T1Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF3T2Y1_prop_f<-asv.n0_ITSF3T2Y1[apply(asv.n0_ITSF3T2Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF3T3Y1_prop_f<-asv.n0_ITSF3T3Y1[apply(asv.n0_ITSF3T3Y1_prop, 1, min) > 0.000001, ]
asv.n0_ITSF3T4Y1_prop_f<-asv.n0_ITSF3T4Y1[apply(asv.n0_ITSF3T4Y1_prop, 1, min) > 0.000001, ]

asv.n0_ITSF1T1Y2_prop_f<-asv.n0_ITSF1T1Y2[apply(asv.n0_ITSF1T1Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF1T2Y2_prop_f<-asv.n0_ITSF1T2Y2[apply(asv.n0_ITSF1T2Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF1T3Y2_prop_f<-asv.n0_ITSF1T3Y2[apply(asv.n0_ITSF1T3Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF1T4Y2_prop_f<-asv.n0_ITSF1T4Y2[apply(asv.n0_ITSF1T4Y2_prop, 1, min) > 0.000001, ]

asv.n0_ITSF2T1Y2_prop_f<-asv.n0_ITSF2T1Y2[apply(asv.n0_ITSF2T1Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF2T2Y2_prop_f<-asv.n0_ITSF2T2Y2[apply(asv.n0_ITSF2T2Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF2T3Y2_prop_f<-asv.n0_ITSF2T3Y2[apply(asv.n0_ITSF2T3Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF2T4Y2_prop_f<-asv.n0_ITSF2T4Y2[apply(asv.n0_ITSF2T4Y2_prop, 1, min) > 0.000001, ]

asv.n0_ITSF3T1Y2_prop_f<-asv.n0_ITSF3T1Y2[apply(asv.n0_ITSF3T1Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF3T2Y2_prop_f<-asv.n0_ITSF3T2Y2[apply(asv.n0_ITSF3T2Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF3T3Y2_prop_f<-asv.n0_ITSF3T3Y2[apply(asv.n0_ITSF3T3Y2_prop, 1, min) > 0.000001, ]
asv.n0_ITSF3T4Y2_prop_f<-asv.n0_ITSF3T4Y2[apply(asv.n0_ITSF3T4Y2_prop, 1, min) > 0.000001, ]

#Check that samples are on columns and ASVs in rows
head(asv.n0_16sF1T1Y1_prop_f) 
head(asv.n0_16sF1T2Y1_prop_f) 
head(asv.n0_16sF1T3Y1_prop_f) 
head(asv.n0_16sF1T4Y1_prop_f) 

head(asv.n0_16sF2T1Y1_prop_f) 
head(asv.n0_16sF2T2Y1_prop_f) 
head(asv.n0_16sF2T3Y1_prop_f) 
head(asv.n0_16sF2T4Y1_prop_f) 

head(asv.n0_16sF3T1Y1_prop_f) 
head(asv.n0_16sF3T2Y1_prop_f) 
head(asv.n0_16sF3T3Y1_prop_f) 
head(asv.n0_16sF3T4Y1_prop_f) 

head(asv.n0_16sF1T1Y2_prop_f) 
head(asv.n0_16sF1T2Y2_prop_f) 
head(asv.n0_16sF1T3Y2_prop_f) 
head(asv.n0_16sF1T4Y2_prop_f) 

head(asv.n0_16sF2T1Y2_prop_f) 
head(asv.n0_16sF2T2Y2_prop_f) 
head(asv.n0_16sF2T3Y2_prop_f) 
head(asv.n0_16sF2T4Y2_prop_f) 

head(asv.n0_16sF3T1Y2_prop_f) 
head(asv.n0_16sF3T2Y2_prop_f) 
head(asv.n0_16sF3T3Y2_prop_f) 
head(asv.n0_16sF3T4Y2_prop_f) 

head(asv.n0_ITSF1T1Y1_prop_f) 
head(asv.n0_ITSF1T2Y1_prop_f) 
head(asv.n0_ITSF1T3Y1_prop_f) 
head(asv.n0_ITSF1T4Y1_prop_f) 

head(asv.n0_ITSF2T1Y1_prop_f) 
head(asv.n0_ITSF2T2Y1_prop_f) 
head(asv.n0_ITSF2T3Y1_prop_f) 
head(asv.n0_ITSF2T4Y1_prop_f) 

head(asv.n0_ITSF3T1Y1_prop_f) 
head(asv.n0_ITSF3T2Y1_prop_f) 
head(asv.n0_ITSF3T3Y1_prop_f) 
head(asv.n0_ITSF3T4Y1_prop_f) 

head(asv.n0_ITSF1T1Y2_prop_f) 
head(asv.n0_ITSF1T2Y2_prop_f) 
head(asv.n0_ITSF1T3Y2_prop_f) 
head(asv.n0_ITSF1T4Y2_prop_f) 

head(asv.n0_ITSF2T1Y2_prop_f) 
head(asv.n0_ITSF2T2Y2_prop_f) 
head(asv.n0_ITSF2T3Y2_prop_f) 
head(asv.n0_ITSF2T4Y2_prop_f) 

head(asv.n0_ITSF3T1Y2_prop_f) 
head(asv.n0_ITSF3T2Y2_prop_f) 
head(asv.n0_ITSF3T3Y2_prop_f) 
head(asv.n0_ITSF3T4Y2_prop_f) 

#Step 5: perform CLR transformation
asv.n0.clr_16sF1T1Y1<-t(apply(asv.n0_16sF1T1Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF1T2Y1<-t(apply(asv.n0_16sF1T2Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF1T3Y1<-t(apply(asv.n0_16sF1T3Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF1T4Y1<-t(apply(asv.n0_16sF1T4Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_16sF2T1Y1<-t(apply(asv.n0_16sF2T1Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF2T2Y1<-t(apply(asv.n0_16sF2T2Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF2T3Y1<-t(apply(asv.n0_16sF2T3Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF2T4Y1<-t(apply(asv.n0_16sF2T4Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_16sF3T1Y1<-t(apply(asv.n0_16sF3T1Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF3T2Y1<-t(apply(asv.n0_16sF3T2Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF3T3Y1<-t(apply(asv.n0_16sF3T3Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF3T4Y1<-t(apply(asv.n0_16sF3T4Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_16sF1T1Y2<-t(apply(asv.n0_16sF1T1Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF1T2Y2<-t(apply(asv.n0_16sF1T2Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF1T3Y2<-t(apply(asv.n0_16sF1T3Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF1T4Y2<-t(apply(asv.n0_16sF1T4Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_16sF2T1Y2<-t(apply(asv.n0_16sF2T1Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF2T2Y2<-t(apply(asv.n0_16sF2T2Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF2T3Y2<-t(apply(asv.n0_16sF2T3Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF2T4Y2<-t(apply(asv.n0_16sF2T4Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_16sF3T1Y2<-t(apply(asv.n0_16sF3T1Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF3T2Y2<-t(apply(asv.n0_16sF3T2Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF3T3Y2<-t(apply(asv.n0_16sF3T3Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_16sF3T4Y2<-t(apply(asv.n0_16sF3T4Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))


asv.n0.clr_ITSF1T1Y1<-t(apply(asv.n0_ITSF1T1Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF1T2Y1<-t(apply(asv.n0_ITSF1T2Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF1T3Y1<-t(apply(asv.n0_ITSF1T3Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF1T4Y1<-t(apply(asv.n0_ITSF1T4Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_ITSF2T1Y1<-t(apply(asv.n0_ITSF2T1Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF2T2Y1<-t(apply(asv.n0_ITSF2T2Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF2T3Y1<-t(apply(asv.n0_ITSF2T3Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF2T4Y1<-t(apply(asv.n0_ITSF2T4Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_ITSF3T1Y1<-t(apply(asv.n0_ITSF3T1Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF3T2Y1<-t(apply(asv.n0_ITSF3T2Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF3T3Y1<-t(apply(asv.n0_ITSF3T3Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF3T4Y1<-t(apply(asv.n0_ITSF3T4Y1_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_ITSF1T1Y2<-t(apply(asv.n0_ITSF1T1Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF1T2Y2<-t(apply(asv.n0_ITSF1T2Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF1T3Y2<-t(apply(asv.n0_ITSF1T3Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF1T4Y2<-t(apply(asv.n0_ITSF1T4Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_ITSF2T1Y2<-t(apply(asv.n0_ITSF2T1Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF2T2Y2<-t(apply(asv.n0_ITSF2T2Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF2T3Y2<-t(apply(asv.n0_ITSF2T3Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF2T4Y2<-t(apply(asv.n0_ITSF2T4Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))

asv.n0.clr_ITSF3T1Y2<-t(apply(asv.n0_ITSF3T1Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF3T2Y2<-t(apply(asv.n0_ITSF3T2Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF3T3Y2<-t(apply(asv.n0_ITSF3T3Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))
asv.n0.clr_ITSF3T4Y2<-t(apply(asv.n0_ITSF3T4Y2_prop_f, 2, function(x){log(x)-mean(log(x))}))

#Check output table. Samples should be in rows and ASVs in columns
head(asv.n0.clr_16sF1T1Y1) 
head(asv.n0.clr_16sF1T2Y1) 
head(asv.n0.clr_16sF1T3Y1) 
head(asv.n0.clr_16sF1T4Y1) 

head(asv.n0.clr_16sF2T1Y1) 
head(asv.n0.clr_16sF2T2Y1) 
head(asv.n0.clr_16sF2T3Y1) 
head(asv.n0.clr_16sF2T4Y1) 

head(asv.n0.clr_16sF3T1Y1) 
head(asv.n0.clr_16sF3T2Y1) 
head(asv.n0.clr_16sF3T3Y1) 
head(asv.n0.clr_16sF3T4Y1) 

head(asv.n0.clr_16sF1T1Y2) 
head(asv.n0.clr_16sF1T2Y2) 
head(asv.n0.clr_16sF1T3Y2) 
head(asv.n0.clr_16sF1T4Y2) 

head(asv.n0.clr_16sF2T1Y2) 
head(asv.n0.clr_16sF2T2Y2) 
head(asv.n0.clr_16sF2T3Y2) 
head(asv.n0.clr_16sF2T4Y2) 

head(asv.n0.clr_16sF3T1Y2) 
head(asv.n0.clr_16sF3T2Y2) 
head(asv.n0.clr_16sF3T3Y2) 
head(asv.n0.clr_16sF3T4Y2) 


head(asv.n0.clr_ITSF1T1Y1) 
head(asv.n0.clr_ITSF1T2Y1) 
head(asv.n0.clr_ITSF1T3Y1) 
head(asv.n0.clr_ITSF1T4Y1) 

head(asv.n0.clr_ITSF2T1Y1) 
head(asv.n0.clr_ITSF2T2Y1) 
head(asv.n0.clr_ITSF2T3Y1) 
head(asv.n0.clr_ITSF2T4Y1) 

head(asv.n0.clr_ITSF3T1Y1) 
head(asv.n0.clr_ITSF3T2Y1) 
head(asv.n0.clr_ITSF3T3Y1) 
head(asv.n0.clr_ITSF3T4Y1) 

head(asv.n0.clr_ITSF1T1Y2) 
head(asv.n0.clr_ITSF1T2Y2) 
head(asv.n0.clr_ITSF1T3Y2) 
head(asv.n0.clr_ITSF1T4Y2) 

head(asv.n0.clr_ITSF2T1Y2) 
head(asv.n0.clr_ITSF2T2Y2) 
head(asv.n0.clr_ITSF2T3Y2) 
head(asv.n0.clr_ITSF2T4Y2) 

head(asv.n0.clr_ITSF3T1Y2) 
head(asv.n0.clr_ITSF3T2Y2) 
head(asv.n0.clr_ITSF3T3Y2) 
head(asv.n0.clr_ITSF3T4Y2) 

#PERMANOVAs
#Calculate Aitchinson distance
dist_16sF1T1Y1<-dist(asv.n0.clr_16sF1T1Y1, method='euclidean')
dist_16sF1T2Y1<-dist(asv.n0.clr_16sF1T2Y1, method='euclidean')
dist_16sF1T3Y1<-dist(asv.n0.clr_16sF1T3Y1, method='euclidean')
dist_16sF1T4Y1<-dist(asv.n0.clr_16sF1T4Y1, method='euclidean')

dist_16sF2T1Y1<-dist(asv.n0.clr_16sF2T1Y1, method='euclidean')
dist_16sF2T2Y1<-dist(asv.n0.clr_16sF2T2Y1, method='euclidean')
dist_16sF2T3Y1<-dist(asv.n0.clr_16sF2T3Y1, method='euclidean')
dist_16sF2T4Y1<-dist(asv.n0.clr_16sF2T4Y1, method='euclidean')

dist_16sF3T1Y1<-dist(asv.n0.clr_16sF3T1Y1, method='euclidean')
dist_16sF3T2Y1<-dist(asv.n0.clr_16sF3T2Y1, method='euclidean')
dist_16sF3T3Y1<-dist(asv.n0.clr_16sF3T3Y1, method='euclidean')
dist_16sF3T4Y1<-dist(asv.n0.clr_16sF3T4Y1, method='euclidean')

dist_16sF1T1Y2<-dist(asv.n0.clr_16sF1T1Y2, method='euclidean')
dist_16sF1T2Y2<-dist(asv.n0.clr_16sF1T2Y2, method='euclidean')
dist_16sF1T3Y2<-dist(asv.n0.clr_16sF1T3Y2, method='euclidean')
dist_16sF1T4Y2<-dist(asv.n0.clr_16sF1T4Y2, method='euclidean')

dist_16sF2T1Y2<-dist(asv.n0.clr_16sF2T1Y2, method='euclidean')
dist_16sF2T2Y2<-dist(asv.n0.clr_16sF2T2Y2, method='euclidean')
dist_16sF2T3Y2<-dist(asv.n0.clr_16sF2T3Y2, method='euclidean')
dist_16sF2T4Y2<-dist(asv.n0.clr_16sF2T4Y2, method='euclidean')

dist_16sF3T1Y2<-dist(asv.n0.clr_16sF3T1Y2, method='euclidean')
dist_16sF3T2Y2<-dist(asv.n0.clr_16sF3T2Y2, method='euclidean')
dist_16sF3T3Y2<-dist(asv.n0.clr_16sF3T3Y2, method='euclidean')
dist_16sF3T4Y2<-dist(asv.n0.clr_16sF3T4Y2, method='euclidean')

dist_ITSF1T1Y1<-dist(asv.n0.clr_ITSF1T1Y1, method='euclidean')
dist_ITSF1T2Y1<-dist(asv.n0.clr_ITSF1T2Y1, method='euclidean')
dist_ITSF1T3Y1<-dist(asv.n0.clr_ITSF1T3Y1, method='euclidean')
dist_ITSF1T4Y1<-dist(asv.n0.clr_ITSF1T4Y1, method='euclidean')

dist_ITSF2T1Y1<-dist(asv.n0.clr_ITSF2T1Y1, method='euclidean')
dist_ITSF2T2Y1<-dist(asv.n0.clr_ITSF2T2Y1, method='euclidean')
dist_ITSF2T3Y1<-dist(asv.n0.clr_ITSF2T3Y1, method='euclidean')
dist_ITSF2T4Y1<-dist(asv.n0.clr_ITSF2T4Y1, method='euclidean')

dist_ITSF3T1Y1<-dist(asv.n0.clr_ITSF3T1Y1, method='euclidean')
dist_ITSF3T2Y1<-dist(asv.n0.clr_ITSF3T2Y1, method='euclidean')
dist_ITSF3T3Y1<-dist(asv.n0.clr_ITSF3T3Y1, method='euclidean')
dist_ITSF3T4Y1<-dist(asv.n0.clr_ITSF3T4Y1, method='euclidean')

dist_ITSF1T1Y2<-dist(asv.n0.clr_ITSF1T1Y2, method='euclidean')
dist_ITSF1T2Y2<-dist(asv.n0.clr_ITSF1T2Y2, method='euclidean')
dist_ITSF1T3Y2<-dist(asv.n0.clr_ITSF1T3Y2, method='euclidean')
dist_ITSF1T4Y2<-dist(asv.n0.clr_ITSF1T4Y2, method='euclidean')

dist_ITSF2T1Y2<-dist(asv.n0.clr_ITSF2T1Y2, method='euclidean')
dist_ITSF2T2Y2<-dist(asv.n0.clr_ITSF2T2Y2, method='euclidean')
dist_ITSF2T3Y2<-dist(asv.n0.clr_ITSF2T3Y2, method='euclidean')
dist_ITSF2T4Y2<-dist(asv.n0.clr_ITSF2T4Y2, method='euclidean')

dist_ITSF3T1Y2<-dist(asv.n0.clr_ITSF3T1Y2, method='euclidean')
dist_ITSF3T2Y2<-dist(asv.n0.clr_ITSF3T2Y2, method='euclidean')
dist_ITSF3T3Y2<-dist(asv.n0.clr_ITSF3T3Y2, method='euclidean')
dist_ITSF3T4Y2<-dist(asv.n0.clr_ITSF3T4Y2, method='euclidean')


#One-way PERMANOVA by treatment_time interaction
permanova_16sF1T1Y1<-pairwise.adonis(dist_16sF1T1Y1, factors = metadata.16sF1T1Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1T2Y1<-pairwise.adonis(dist_16sF1T2Y1, factors = metadata.16sF1T2Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1T3Y1<-pairwise.adonis(dist_16sF1T3Y1, factors = metadata.16sF1T3Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1T4Y1<-pairwise.adonis(dist_16sF1T4Y1, factors = metadata.16sF1T4Y1$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_16sF2T1Y1<-pairwise.adonis(dist_16sF2T1Y1, factors = metadata.16sF2T1Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2T2Y1<-pairwise.adonis(dist_16sF2T2Y1, factors = metadata.16sF2T2Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2T3Y1<-pairwise.adonis(dist_16sF2T3Y1, factors = metadata.16sF2T3Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2T4Y1<-pairwise.adonis(dist_16sF2T4Y1, factors = metadata.16sF2T4Y1$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_16sF3T1Y1<-pairwise.adonis(dist_16sF3T1Y1, factors = metadata.16sF3T1Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3T2Y1<-pairwise.adonis(dist_16sF3T2Y1, factors = metadata.16sF3T2Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3T3Y1<-pairwise.adonis(dist_16sF3T3Y1, factors = metadata.16sF3T3Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3T4Y1<-pairwise.adonis(dist_16sF3T4Y1, factors = metadata.16sF3T4Y1$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_16sF1T1Y2<-pairwise.adonis(dist_16sF1T1Y2, factors = metadata.16sF1T1Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1T2Y2<-pairwise.adonis(dist_16sF1T2Y2, factors = metadata.16sF1T2Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1T3Y2<-pairwise.adonis(dist_16sF1T3Y2, factors = metadata.16sF1T3Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1T4Y2<-pairwise.adonis(dist_16sF1T4Y2, factors = metadata.16sF1T4Y2$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_16sF2T1Y2<-pairwise.adonis(dist_16sF2T1Y2, factors = metadata.16sF2T1Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2T2Y2<-pairwise.adonis(dist_16sF2T2Y2, factors = metadata.16sF2T2Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2T3Y2<-pairwise.adonis(dist_16sF2T3Y2, factors = metadata.16sF2T3Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2T4Y2<-pairwise.adonis(dist_16sF2T4Y2, factors = metadata.16sF2T4Y2$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_16sF3T1Y2<-pairwise.adonis(dist_16sF3T1Y2, factors = metadata.16sF3T1Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3T2Y2<-pairwise.adonis(dist_16sF3T2Y2, factors = metadata.16sF3T2Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3T3Y2<-pairwise.adonis(dist_16sF3T3Y2, factors = metadata.16sF3T3Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3T4Y2<-pairwise.adonis(dist_16sF3T4Y2, factors = metadata.16sF3T4Y2$Time, perm = 999, p.adjust.m = 'bonferroni')


permanova_ITSF1T1Y1<-pairwise.adonis(dist_ITSF1T1Y1, factors = metadata.ITSF1T1Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1T2Y1<-pairwise.adonis(dist_ITSF1T2Y1, factors = metadata.ITSF1T2Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1T3Y1<-pairwise.adonis(dist_ITSF1T3Y1, factors = metadata.ITSF1T3Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1T4Y1<-pairwise.adonis(dist_ITSF1T4Y1, factors = metadata.ITSF1T4Y1$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_ITSF2T1Y1<-pairwise.adonis(dist_ITSF2T1Y1, factors = metadata.ITSF2T1Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2T2Y1<-pairwise.adonis(dist_ITSF2T2Y1, factors = metadata.ITSF2T2Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2T3Y1<-pairwise.adonis(dist_ITSF2T3Y1, factors = metadata.ITSF2T3Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2T4Y1<-pairwise.adonis(dist_ITSF2T4Y1, factors = metadata.ITSF2T4Y1$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_ITSF3T1Y1<-pairwise.adonis(dist_ITSF3T1Y1, factors = metadata.ITSF3T1Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3T2Y1<-pairwise.adonis(dist_ITSF3T2Y1, factors = metadata.ITSF3T2Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3T3Y1<-pairwise.adonis(dist_ITSF3T3Y1, factors = metadata.ITSF3T3Y1$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3T4Y1<-pairwise.adonis(dist_ITSF3T4Y1, factors = metadata.ITSF3T4Y1$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_ITSF1T1Y2<-pairwise.adonis(dist_ITSF1T1Y2, factors = metadata.ITSF1T1Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1T2Y2<-pairwise.adonis(dist_ITSF1T2Y2, factors = metadata.ITSF1T2Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1T3Y2<-pairwise.adonis(dist_ITSF1T3Y2, factors = metadata.ITSF1T3Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1T4Y2<-pairwise.adonis(dist_ITSF1T4Y2, factors = metadata.ITSF1T4Y2$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_ITSF2T1Y2<-pairwise.adonis(dist_ITSF2T1Y2, factors = metadata.ITSF2T1Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2T2Y2<-pairwise.adonis(dist_ITSF2T2Y2, factors = metadata.ITSF2T2Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2T3Y2<-pairwise.adonis(dist_ITSF2T3Y2, factors = metadata.ITSF2T3Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2T4Y2<-pairwise.adonis(dist_ITSF2T4Y2, factors = metadata.ITSF2T4Y2$Time, perm = 999, p.adjust.m = 'bonferroni')

permanova_ITSF3T1Y2<-pairwise.adonis(dist_ITSF3T1Y2, factors = metadata.ITSF3T1Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3T2Y2<-pairwise.adonis(dist_ITSF3T2Y2, factors = metadata.ITSF3T2Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3T3Y2<-pairwise.adonis(dist_ITSF3T3Y2, factors = metadata.ITSF3T3Y2$Time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3T4Y2<-pairwise.adonis(dist_ITSF3T4Y2, factors = metadata.ITSF3T4Y2$Time, perm = 999, p.adjust.m = 'bonferroni')


#Step 6: Perform Singular Value Decomposition (PCA)
pc.clr_16sF1<-prcomp(asv.n0.clr_16sF1)
pc.clr_16sF2<-prcomp(asv.n0.clr_16sF2)
pc.clr_16sF3<-prcomp(asv.n0.clr_16sF3)

pc.clr_ITSF1<-prcomp(asv.n0.clr_ITSF1)
pc.clr_ITSF2<-prcomp(asv.n0.clr_ITSF2)
pc.clr_ITSF3<-prcomp(asv.n0.clr_ITSF3)

png("Screeplot by facility- PCA.png", width = 400, height = 300, units = 'px')
par(mar=c(2,2,2,2))
par(mfrow=c(2,3))
screeplot(pc.clr_16sF1, type='barplot', main="Bacteria F1")
screeplot(pc.clr_16sF2, type='barplot', main="Bacteria F2")
screeplot(pc.clr_16sF3, type='barplot', main="Bacteria F3")
screeplot(pc.clr_ITSF1, type='barplot', main="Fungi F1")
screeplot(pc.clr_ITSF2, type='barplot', main="Fungi F2")
screeplot(pc.clr_ITSF3, type='barplot', main="Fungi F3")
dev.off()

#Calculate total variance of the data
mvar.clr_16sF1<-mvar(asv.n0.clr_16sF1)
mvar.clr_16sF2<-mvar(asv.n0.clr_16sF2)
mvar.clr_16sF3<-mvar(asv.n0.clr_16sF3)

mvar.clr_ITSF1<-mvar(asv.n0.clr_ITSF1)
mvar.clr_ITSF2<-mvar(asv.n0.clr_ITSF2)
mvar.clr_ITSF3<-mvar(asv.n0.clr_ITSF3)

#Display results - 16s
row_16sF1<-rownames(asv.n0.clr_16sF1) #Make vector with sample names
pc_out_16sF1<-as.data.frame(pc.clr_16sF1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF1<-as.data.frame(bind_cols(pc_out_16sF1,metadata.16sF1)) #Add metadata information
row.names(pc_out_meta_16sF1)<-row_16sF1 #Add rownames to dataframe
pc_out_meta_16sF1$Facility<-as.factor(pc_out_meta_16sF1$Facility)
pc_out_meta_16sF1$Time<-as.factor(pc_out_meta_16sF1$Time)
pc_out_meta_16sF1$Treatment<-as.factor(pc_out_meta_16sF1$Treatment)

row_16sF2<-rownames(asv.n0.clr_16sF2) #Make vector with sample names
pc_out_16sF2<-as.data.frame(pc.clr_16sF2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF2<-as.data.frame(bind_cols(pc_out_16sF2,metadata.16sF2)) #Add metadata information
row.names(pc_out_meta_16sF2)<-row_16sF2 #Add rownames to dataframe
pc_out_meta_16sF2$Facility<-as.factor(pc_out_meta_16sF2$Facility)
pc_out_meta_16sF2$Time<-as.factor(pc_out_meta_16sF2$Time)
pc_out_meta_16sF2$Treatment<-as.factor(pc_out_meta_16sF2$Treatment)

row_16sF3<-rownames(asv.n0.clr_16sF3) #Make vector with sample names
pc_out_16sF3<-as.data.frame(pc.clr_16sF3$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_16sF3<-as.data.frame(bind_cols(pc_out_16sF3,metadata.16sF3)) #Add metadata information
row.names(pc_out_meta_16sF3)<-row_16sF3 #Add rownames to dataframe
pc_out_meta_16sF3$Facility<-as.factor(pc_out_meta_16sF3$Facility)
pc_out_meta_16sF3$Time<-as.factor(pc_out_meta_16sF3$Time)
pc_out_meta_16sF3$Treatment<-as.factor(pc_out_meta_16sF3$Treatment)

# Make PCA plot - First 2 axis- color by facility/Y and shape by L..mono
#Fig 2A
PCA_16sF1 <- ggplot(pc_out_meta_16sF1, aes(x=PC1,y=PC2, color=Treatment, shape=Time))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF1$sdev[1]^2/mvar.clr_16sF1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF1$sdev[2]^2/mvar.clr_16sF1*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria F1", subtitle = "PCA by Treatment and Time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(option='viridis')
PCA_16sF1
ggsave("PCA_Bacteria_F1_ASV.png", plot =PCA_16sF1, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Bacteria_F1_ASV.svg", plot =PCA_16sF1, device="svg", width=6, height=5, units="in",dpi=600)

PCA_16sF2 <- ggplot(pc_out_meta_16sF2, aes(x=PC1,y=PC2, color=Treatment, shape=Time))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF2$sdev[1]^2/mvar.clr_16sF2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF2$sdev[2]^2/mvar.clr_16sF2*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria F2", subtitle = "PCA by Treatment and Time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(option='viridis')
PCA_16sF2
ggsave("PCA_Bacteria_F2_ASV.png", plot =PCA_16sF2, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Bacteria_F2_ASV.svg", plot =PCA_16sF2, device="svg", width=6, height=5, units="in",dpi=600)

PCA_16sF3 <- ggplot(pc_out_meta_16sF3, aes(x=PC1,y=PC2, color=Treatment, shape=Time))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_16sF3$sdev[1]^2/mvar.clr_16sF3*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_16sF3$sdev[2]^2/mvar.clr_16sF3*100, digits=1), "%", sep="")) +
  ggtitle("Bacteria F3", subtitle = "PCA by Treatment and Time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(option='viridis')
PCA_16sF3
ggsave("PCA_Bacteria_F3_ASV.png", plot =PCA_16sF3, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Bacteria_F3_ASV.svg", plot =PCA_16sF3, device="svg", width=6, height=5, units="in",dpi=600)


#Display results - ITS
row_ITSF1<-rownames(asv.n0.clr_ITSF1) #Make vector with sample names
pc_out_ITSF1<-as.data.frame(pc.clr_ITSF1$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF1<-as.data.frame(bind_cols(pc_out_ITSF1,metadata.ITSF1)) #Add metadata information
row.names(pc_out_meta_ITSF1)<-row_ITSF1 #Add rownames to dataframe
pc_out_meta_ITSF1$Facility<-as.factor(pc_out_meta_ITSF1$Facility)
pc_out_meta_ITSF1$Time<-as.factor(pc_out_meta_ITSF1$Time)
pc_out_meta_ITSF1$Treatment<-as.factor(pc_out_meta_ITSF1$Treatment)

row_ITSF2<-rownames(asv.n0.clr_ITSF2) #Make vector with sample names
pc_out_ITSF2<-as.data.frame(pc.clr_ITSF2$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF2<-as.data.frame(bind_cols(pc_out_ITSF2,metadata.ITSF2)) #Add metadata information
row.names(pc_out_meta_ITSF2)<-row_ITSF2 #Add rownames to dataframe
pc_out_meta_ITSF2$Facility<-as.factor(pc_out_meta_ITSF2$Facility)
pc_out_meta_ITSF2$Time<-as.factor(pc_out_meta_ITSF2$Time)
pc_out_meta_ITSF2$Treatment<-as.factor(pc_out_meta_ITSF2$Treatment)

row_ITSF3<-rownames(asv.n0.clr_ITSF3) #Make vector with sample names
pc_out_ITSF3<-as.data.frame(pc.clr_ITSF3$x[,1:2]) #Get PC1 and PC2 
pc_out_meta_ITSF3<-as.data.frame(bind_cols(pc_out_ITSF3,metadata.ITSF3)) #Add metadata information
row.names(pc_out_meta_ITSF3)<-row_ITSF3 #Add rownames to dataframe
pc_out_meta_ITSF3$Facility<-as.factor(pc_out_meta_ITSF3$Facility)
pc_out_meta_ITSF3$Time<-as.factor(pc_out_meta_ITSF3$Time)
pc_out_meta_ITSF3$Treatment<-as.factor(pc_out_meta_ITSF3$Treatment)

# Make PCA plot - First 2 axis- color by facility/Y and shape by L..mono
#Fig 2A
PCA_ITSF1 <- ggplot(pc_out_meta_ITSF1, aes(x=PC1,y=PC2, color=Treatment, shape=Time))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF1$sdev[1]^2/mvar.clr_ITSF1*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF1$sdev[2]^2/mvar.clr_ITSF1*100, digits=1), "%", sep="")) +
  ggtitle("Fungi F1", subtitle = "PCA by Treatment and Time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(option='viridis')
PCA_ITSF1
ggsave("PCA_Fungi_F1_ASV.png", plot =PCA_ITSF1, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Fungi_F1_ASV.svg", plot =PCA_ITSF1, device="svg", width=6, height=5, units="in",dpi=600)

PCA_ITSF2 <- ggplot(pc_out_meta_ITSF2, aes(x=PC1,y=PC2, color=Treatment, shape=Time))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF2$sdev[1]^2/mvar.clr_ITSF2*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF2$sdev[2]^2/mvar.clr_ITSF2*100, digits=1), "%", sep="")) +
  ggtitle("Fungi F2", subtitle = "PCA by Treatment and Time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(option='viridis')
PCA_ITSF2
ggsave("PCA_Fungi_F2_ASV.png", plot =PCA_ITSF2, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Fungi_F2_ASV.svg", plot =PCA_ITSF2, device="svg", width=6, height=5, units="in",dpi=600)

PCA_ITSF3 <- ggplot(pc_out_meta_ITSF3, aes(x=PC1,y=PC2, color=Treatment, shape=Time))+
  geom_point(size=3)+
  theme(legend.text=element_text(size=15,color='black'), legend.title= element_text(size=15, face='italic')) +
  theme(legend.position = 'right')+
  #theme(plot.margin=margin(t=0.5, b=0.5, l=0.5, r=0.5, unit = 'in'))+
  theme(axis.text = element_blank(), axis.ticks = element_blank()) +
  theme(axis.title = element_text(size=15,color='black')) +
  theme(panel.background = element_rect(fill='white', color = NA),
        plot.background = element_rect(fill = 'white',color = NA), 
        panel.border = element_rect(color='black',fill = NA,size=1))+
  scale_x_continuous(name = paste("PC1: ", round(pc.clr_ITSF3$sdev[1]^2/mvar.clr_ITSF3*100, digits=1), "%", sep="")) +
  scale_y_continuous(name = paste("PC2: ", round(pc.clr_ITSF3$sdev[2]^2/mvar.clr_ITSF3*100, digits=1), "%", sep="")) +
  ggtitle("Fungi F3", subtitle = "PCA by Treatment and Time")+theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=16))+
  scale_color_viridis_d(option='viridis')
PCA_ITSF3
ggsave("PCA_Fungi_F3_ASV.png", plot =PCA_ITSF3, device="png", width=6, height=5, units="in",dpi=600)
ggsave("PCA_Fungi_F3_ASV.svg", plot =PCA_ITSF3, device="svg", width=6, height=5, units="in",dpi=600)

# PERMANOVA #
#Calculate Aitchinson distance
dist_16sF1<-dist(t(asv.n0.clr_16sF1), method='euclidean')
dist_16sF2<-dist(t(asv.n0.clr_16sF2), method='euclidean')
dist_16sF3<-dist(t(asv.n0.clr_16sF3), method='euclidean')

dist_ITSF1<-dist(t(asv.n0.clr_ITSF1), method='euclidean')
dist_ITSF2<-dist(t(asv.n0.clr_ITSF2), method='euclidean')
dist_ITSF3<-dist(t(asv.n0.clr_ITSF3), method='euclidean')


#One-way PERMANOVA by treatment_time interaction
permanova_16sF1_trt<-pairwise.adonis(dist_16sF1, factors = metadata.16sF1$Trt_time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1_trt

permanova_16sF2_trt<-pairwise.adonis(dist_16sF2, factors = metadata.16sF2$Trt_time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2_trt

permanova_16sF3_trt<-pairwise.adonis(dist_16sF3, factors = metadata.16sF3$Trt_time, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3_trt

permanova_ITSF1_trt<-pairwise.adonis(dist_ITSF1, factors = metadata.ITSF1$Trt_time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1_trt

permanova_ITSF2_trt<-pairwise.adonis(dist_ITSF2, factors = metadata.ITSF2$Trt_time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2_trt

permanova_ITSF3_trt<-pairwise.adonis(dist_ITSF3, factors = metadata.ITSF3$Trt_time, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3_trt


#Two way ANOVA by Treatment and Time
permanova_16sF1<-pairwise.adonis2(dist_16sF1~Treatment+Time, data=metadata.16sF1, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF1

permanova_16sF2<-pairwise.adonis2(dist_16sF2~Treatment+Time, data=metadata.16sF2, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF2

permanova_16sF3<-pairwise.adonis2(dist_16sF3~Treatment+Time, data=metadata.16sF3, perm = 999, p.adjust.m = 'bonferroni')
permanova_16sF3

#ITS
permanova_ITSF1<-pairwise.adonis2(dist_ITSF1~+Treatment+Time, data=metadata.ITSF1, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF1

permanova_ITSF2<-pairwise.adonis2(dist_ITSF2~+Treatment+Time, data=metadata.ITSF2, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF2

permanova_ITSF3<-pairwise.adonis2(dist_ITSF3~Treatment+Time, data=metadata.ITSF3, perm = 999, p.adjust.m = 'bonferroni')
permanova_ITSF3


##NOTE: Permanova showed significant difference by treatment for all bacteria datasets. for fungi, there were sig difference by teatment in F1
#by trt and time in F2 and by trrt, time and interaction in F3


# 
# 
# 
# 
# 
# 
# #### Differential Abundance analysis by Facility ####
# library(ALDEx2)
# 
# #Change Time to character vector - as.factor affects aldex.effect function
# metadata.16sF1$Time <- as.character(metadata.16sF1$Time)
# metadata.16sF2$Time <- as.character(metadata.16sF2$Time)
# metadata.16sF3$Time <- as.character(metadata.16sF3$Time)
# 
# metadata.ITSF1$Time <- as.character(metadata.ITSF1$Time)
# metadata.ITSF2$Time <- as.character(metadata.ITSF2$Time)
# metadata.ITSF3$Time <- as.character(metadata.ITSF3$Time)
# 
# #ASV table needs to be as ASV in rows. O
# asv_16sF1<-t(asv_16sF1)
# asv_16sF2<-t(asv_16sF2)
# asv_16sF3<-t(asv_16sF3)
# 
# asv_ITSF1<-t(asv_ITSF1)
# asv_ITSF2<-t(asv_ITSF2)
# asv_ITSF3<-t(asv_ITSF3)
# 
# #Aldex2
# #Generate 128 Dirichlet distributed Monte Carlo instances and center-log ratio transform them. nly two conditions can be compared at a time.
# Aldex_16sF1.clr<-aldex.clr(asv_16sF1, mc.samples = 128, conds = metadata.16sF1$Time)
# Aldex_16sF2.clr<-aldex.clr(asv_16sF2, mc.samples = 128, conds = metadata.16sF2$Time)
# Aldex_16sF3.clr<-aldex.clr(asv_16sF3, mc.samples = 128, conds = metadata.16sF3$Time)
# 
# Aldex_ITSF1.clr<-aldex.clr(asv_ITSF1, mc.samples = 128, conds = metadata.ITSF1$Time)
# Aldex_ITSF2.clr<-aldex.clr(asv_ITSF2, mc.samples = 128, conds = metadata.ITSF2$Time)
# Aldex_ITSF3.clr<-aldex.clr(asv_ITSF3, mc.samples = 128, conds = metadata.ITSF3$Time)
# 
# #Calculate the expected effect size
# Aldex_16sF1.e<-aldex.effect(Aldex_16sF1.clr)
# Aldex_16sF2.e<-aldex.effect(Aldex_16sF2.clr)
# Aldex_16sF3.e<-aldex.effect(Aldex_16sF3.clr)
# 
# Aldex_ITSF1.e<-aldex.effect(Aldex_ITSF1.clr)
# Aldex_ITSF2.e<-aldex.effect(Aldex_ITSF2.clr)
# Aldex_ITSF3.e<-aldex.effect(Aldex_ITSF3.clr)
# 
# #Generate the P and Benjamini-Hochberg corrected P values.
# Aldex_16sF1.t<-aldex.ttest(Aldex_16sF1.clr)
# Aldex_16sF2.t<-aldex.ttest(Aldex_16sF2.clr)
# Aldex_16sF3.t<-aldex.ttest(Aldex_16sF3.clr)
# 
# Aldex_ITSF1.t<-aldex.ttest(Aldex_ITSF1.clr)
# Aldex_ITSF2.t<-aldex.ttest(Aldex_ITSF2.clr)
# Aldex_ITSF3.t<-aldex.ttest(Aldex_ITSF3.clr)
# 
# #Merge data frames
# Aldex_16sF1.all<-data.frame(Aldex_16sF1.e,Aldex_16sF1.t)
# Aldex_16sF2.all<-data.frame(Aldex_16sF2.e,Aldex_16sF2.t)
# Aldex_16sF3.all<-data.frame(Aldex_16sF3.e,Aldex_16sF3.t)
# 
# Aldex_ITSF1.all<-data.frame(Aldex_ITSF1.e,Aldex_ITSF1.t)
# Aldex_ITSF2.all<-data.frame(Aldex_ITSF2.e,Aldex_ITSF2.t)
# Aldex_ITSF3.all<-data.frame(Aldex_ITSF3.e,Aldex_ITSF3.t)
# 
# #Determine which corrected values fall below a threshold
# Aldex_16sF1.sig<-which(Aldex_16sF1.all$wi.eBH <=0.05)
# Aldex_16sF2.sig<-which(Aldex_16sF2.all$wi.eBH <=0.05)
# Aldex_16sF3.sig<-which(Aldex_16sF3.all$wi.eBH <=0.05)
# 
# Aldex_ITSF1.sig<-which(Aldex_ITSF1.all$wi.eBH <=0.05)
# Aldex_ITSF2.sig<-which(Aldex_ITSF2.all$wi.eBH <=0.05)
# Aldex_ITSF3.sig<-which(Aldex_ITSF3.all$wi.eBH <=0.05)
# 
# 
# #Plots of significant families relative abundances
# #Extract significant ASV data from ALDEx2 output
# # 16s
# Aldex_ITSF2.sig.row<-rownames(Aldex_ITSF2.all)[which(Aldex_ITSF2.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_ITSF2.sig.table<-subset(Aldex_ITSF2.all, rownames(Aldex_ITSF2.all) %in% Aldex_ITSF2.sig.row) #Subset significant families
# Aldex_ITSF2.sig.taxon<-subset(taxon.ITS_clean, rownames(taxon.ITS_clean) %in% Aldex_ITSF2.sig.row) #Subset taxon table to get the taxonomy of significant families
# Aldex_ITSF2.sig.table.all<-bind_cols(Aldex_ITSF2.sig.taxon, Aldex_ITSF2.sig.table) #combine tables
# Aldex_ITSF2.sig.table.all$OTU<-rownames(Aldex_ITSF2.sig.table.all)
# 
# Aldex_ITSF3.sig.row<-rownames(Aldex_ITSF3.all)[which(Aldex_ITSF3.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_ITSF3.sig.table<-subset(Aldex_ITSF3.all, rownames(Aldex_ITSF3.all) %in% Aldex_ITSF3.sig.row) #Subset significant families
# Aldex_ITSF3.sig.taxon<-subset(taxon.ITS_clean, rownames(taxon.ITS_clean) %in% Aldex_ITSF3.sig.row) #Subset taxon table to get the taxonomy of significant families
# Aldex_ITSF3.sig.table.all<-bind_cols(Aldex_ITSF3.sig.taxon, Aldex_ITSF3.sig.table) #combine tables
# Aldex_ITSF3.sig.table.all$OTU<-rownames(Aldex_ITSF3.sig.table.all)
# 
#  
# 
# #Effect plots for those ASV with effect size over 1 and less than -1
# Aldex_ITSF2.eff.table.all<-bind_rows(subset(Aldex_ITSF2.sig.table.all, effect >=1),subset(Aldex_ITSF2.sig.table.all, effect <=-1))
# Aldex_ITSF3.eff.table.all<-bind_rows(subset(Aldex_ITSF3.sig.table.all, effect >=1),subset(Aldex_ITSF3.sig.table.all, effect <=-1))
# 
# 
# ## Boxplots of ASVs with effect size grater than 1
# 
# asv.n0.acomp_ITSF2<-as.data.frame(acomp(t(asv_n0_ITSF2)), total=1)
# asv.n0.acomp_ITSF3<-as.data.frame(acomp(t(asv_n0_ITSF3)), total=1)
# 
# #ASV level plots
# 
# #Make Phyloseq object
# phyloseqITSF2 <- phyloseq(otu_table(asv.n0.acomp_ITSF2, taxa_are_rows = TRUE), tax_table(taxon_ITS), sample_data(metadata.ITSF2))
# phyloseqITSF3 <- phyloseq(otu_table(asv.n0.acomp_ITSF3, taxa_are_rows = TRUE), tax_table(taxon_ITS), sample_data(metadata.ITSF3))
# 
# 
# #ASV level Plots
# #Make long format table from Phyloseq object
# asv_ITSF2_long <- phyloseqITSF2 %>%  
#   transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
#   psmelt() %>%  #Melts to long format
#   arrange(desc(Abundance))
# 
# asv_ITSF3_long <- phyloseqITSF3 %>%  
#   transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
#   psmelt() %>%  #Melts to long format
#   arrange(desc(Abundance))
# 
# 
# #Make vector with ASV names for differential abundant ASVs with effect size over 1 for each pair of comparisons
# asv_ITSF2_effect<-Aldex_ITSF2.eff.table.all$OTU
# asv_ITSF3_effect<-Aldex_ITSF3.eff.table.all$OTU
# 
# 
# #Filter table to obtain only ASVs that with effect size greater than absolute 1
# asv_ITSF2_filter_effect <- filter(asv_ITSF2_long, OTU %in% asv_ITSF2_effect)
# asv_ITSF3_filter_effect <- filter(asv_ITSF3_long, OTU %in% asv_ITSF3_effect)
# 
# 
# #calculate statistics by ASV and Time
# stat_ITSF2<-asv_ITSF2_filter_effect %>%
#   group_by(OTU, Time, Genus)%>%
#   summarise(Mean=mean(Abundance, na.rm=TRUE))
# 
# stat_ITSF3<-asv_ITSF3_filter_effect %>%
#   group_by(OTU, Time, Genus)%>%
#   summarise(Mean=mean(Abundance, na.rm=TRUE))
# 
# 
# #Reshape table to wide format
# stat_ITSF2_mean<-dcast(stat_ITSF2, OTU + Genus~ Time , value.var='Mean')
# stat_ITSF3_mean<-dcast(stat_ITSF3, OTU + Genus~ Time , value.var='Mean')
# 
# #Calculate log fold chage - log2 of the ratio of mean RA of first facility by second facility in comparison
# stat_ITSF2_mean$logchange<-log2(stat_ITSF2_mean$Before/stat_ITSF2_mean$After)
# stat_ITSF3_mean$logchange<-log2(stat_ITSF3_mean$Before/stat_ITSF3_mean$After)
# 
# #Add facility identifier - Year with higher RA 
# stat_ITSF2_mean$Time<-ifelse(stat_ITSF2_mean$logchange <0 , "After", "Before")
# stat_ITSF3_mean$Time<-ifelse(stat_ITSF3_mean$logchange <0 , "After", "Before")
# 
# 
# #Plot 
# Logfold_ITSF2<-ggplot(stat_ITSF2_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Time))+
#   geom_bar(stat='identity', color='black')+
#   geom_text(aes(label=Genus), position = position_stack(vjust=0.5), size=5)+
#   theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=8, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   xlab("Log fold change (log2 Before/After") + 
#   scale_x_continuous(limits=c(-12,12), breaks = c(-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12))+
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="white", color =NA),
#         panel.border = element_rect(color="black", fill=NA)) +
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_manual(values=c("#BB3754","#D6879880"))+
#   #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   ggtitle("Fungi F2 - Before v After", subtitle = "Differential abundant ASVs")
# 
# Logfold_ITSF3<-ggplot(stat_ITSF3_mean, aes(x=logchange, y=reorder(OTU,logchange), fill=Time))+
#   geom_bar(stat='identity', color='black')+
#   geom_text(aes(label=Genus), position = position_stack(vjust=0.5), size=5)+
#   theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=8, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   xlab("Log fold change (log2 Before/After)") + 
#   scale_x_continuous(limits=c(-12,12), breaks = c(-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12))+
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="white", color =NA),
#         panel.border = element_rect(color="black", fill=NA)) +
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_manual(values=c("#FCA50A","#FDC96C80"))+
#   #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   ggtitle("Fungi F3 - Before v After", subtitle = "Differential abundant ASVs")
# 
# 
# 
# #Combine plots
# Logfold_byTime <- plot_grid(Logfold_ITSF2, Logfold_ITSF3,
#                             ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
# Logfold_byTime
# ggsave("Aldex_LogFold_byTime.png", plot=Logfold_byTime, device="png", width=8, height=5, units="in", dpi=600)
# ggsave("Aldex_LogFold_byTime.svg", plot=Logfold_byTime, device="svg", width=8, height=5, units="in", dpi=600)
# 
# 
# #### Differential abundance analysis by facility at GENUS level ####
# ps_16sF1_genus<-tax_glom(ps.16sF1, taxrank = "Genus")
# ps_16sF2_genus<-tax_glom(ps.16sF2, taxrank = "Genus")
# ps_16sF3_genus<-tax_glom(ps.16sF3, taxrank = "Genus")
# 
# ps_ITSF1_genus<-tax_glom(ps.ITSF1, taxrank = "Genus")
# ps_ITSF2_genus<-tax_glom(ps.ITSF2, taxrank = "Genus")
# ps_ITSF3_genus<-tax_glom(ps.ITSF3, taxrank = "Genus")
# 
# 
# #Take otu table for each facility from Phyloseq
# asv_16sF1_genus<-as.data.frame(as(otu_table(ps_16sF1_genus), "matrix"))
# asv_16sF2_genus<-as.data.frame(as(otu_table(ps_16sF2_genus), "matrix"))
# asv_16sF3_genus<-as.data.frame(as(otu_table(ps_16sF3_genus), "matrix"))
# 
# asv_ITSF1_genus<-as.data.frame(as(otu_table(ps_ITSF1_genus), "matrix"))
# asv_ITSF2_genus<-as.data.frame(as(otu_table(ps_ITSF2_genus), "matrix"))
# asv_ITSF3_genus<-as.data.frame(as(otu_table(ps_ITSF3_genus), "matrix"))
# 
# 
# #Transpose table
# asv_16sF1_genus<-t(asv_16sF1_genus)
# asv_16sF2_genus<-t(asv_16sF2_genus)
# asv_16sF3_genus<-t(asv_16sF3_genus)
# 
# asv_ITSF1_genus<-t(asv_ITSF1_genus)
# asv_ITSF2_genus<-t(asv_ITSF2_genus)
# asv_ITSF3_genus<-t(asv_ITSF3_genus)
# 
# #Remove ASVs that have count zero in all samples - Necessary step for the zero imputation function
# asv_16sF1_genus<-asv_16sF1_genus[ which(rowSums(asv_16sF1_genus)>0),]
# asv_16sF2_genus<-asv_16sF2_genus[ which(rowSums(asv_16sF2_genus)>0),]
# asv_16sF3_genus<-asv_16sF3_genus[ which(rowSums(asv_16sF3_genus)>0),]
# 
# asv_ITSF1_genus<-asv_ITSF1_genus[ which(rowSums(asv_ITSF1_genus)>0),]
# asv_ITSF2_genus<-asv_ITSF2_genus[ which(rowSums(asv_ITSF2_genus)>0),]
# asv_ITSF3_genus<-asv_ITSF3_genus[ which(rowSums(asv_ITSF3_genus)>0),]
# 
# 
# #Calculate CLR
# Aldex_16sF1_genus.clr<-aldex.clr(asv_16sF1_genus, mc.samples = 128, conds = metadata.16sF1$Time)
# Aldex_16sF2_genus.clr<-aldex.clr(asv_16sF2_genus, mc.samples = 128, conds = metadata.16sF2$Time)
# Aldex_16sF3_genus.clr<-aldex.clr(asv_16sF3_genus, mc.samples = 128, conds = metadata.16sF3$Time)
# 
# Aldex_ITSF1_genus.clr<-aldex.clr(asv_ITSF1_genus, mc.samples = 128, conds = metadata.ITSF1$Time)
# Aldex_ITSF2_genus.clr<-aldex.clr(asv_ITSF2_genus, mc.samples = 128, conds = metadata.ITSF2$Time)
# Aldex_ITSF3_genus.clr<-aldex.clr(asv_ITSF3_genus, mc.samples = 128, conds = metadata.ITSF3$Time)
# 
# #Calculate the expected effect size
# Aldex_16sF1_genus.e<-aldex.effect(Aldex_16sF1_genus.clr)
# Aldex_16sF2_genus.e<-aldex.effect(Aldex_16sF2_genus.clr)
# Aldex_16sF3_genus.e<-aldex.effect(Aldex_16sF3_genus.clr)
# 
# Aldex_ITSF1_genus.e<-aldex.effect(Aldex_ITSF1_genus.clr)
# Aldex_ITSF2_genus.e<-aldex.effect(Aldex_ITSF2_genus.clr)
# Aldex_ITSF3_genus.e<-aldex.effect(Aldex_ITSF3_genus.clr)
# 
# #Generate the P and Benjamini-Hochberg corrected P values.
# Aldex_16sF1_genus.t<-aldex.ttest(Aldex_16sF1_genus.clr)
# Aldex_16sF2_genus.t<-aldex.ttest(Aldex_16sF2_genus.clr)
# Aldex_16sF3_genus.t<-aldex.ttest(Aldex_16sF3_genus.clr)
# 
# Aldex_ITSF1_genus.t<-aldex.ttest(Aldex_ITSF1_genus.clr)
# Aldex_ITSF2_genus.t<-aldex.ttest(Aldex_ITSF2_genus.clr)
# Aldex_ITSF3_genus.t<-aldex.ttest(Aldex_ITSF3_genus.clr)
# 
# #Merge data frames
# Aldex_16sF1_genus.all<-data.frame(Aldex_16sF1_genus.e,Aldex_16sF1_genus.t)
# Aldex_16sF2_genus.all<-data.frame(Aldex_16sF2_genus.e,Aldex_16sF2_genus.t)
# Aldex_16sF3_genus.all<-data.frame(Aldex_16sF3_genus.e,Aldex_16sF3_genus.t)
# 
# Aldex_ITSF1_genus.all<-data.frame(Aldex_ITSF1_genus.e,Aldex_ITSF1_genus.t)
# Aldex_ITSF2_genus.all<-data.frame(Aldex_ITSF2_genus.e,Aldex_ITSF2_genus.t)
# Aldex_ITSF3_genus.all<-data.frame(Aldex_ITSF3_genus.e,Aldex_ITSF3_genus.t)
# 
# #Determine which corrected values fall below a threshold
# Aldex_16sF1_genus.sig<-which(Aldex_16sF1_genus.all$wi.eBH <=0.05)
# Aldex_16sF2_genus.sig<-which(Aldex_16sF2_genus.all$wi.eBH <=0.05)
# Aldex_16sF3_genus.sig<-which(Aldex_16sF3_genus.all$wi.eBH <=0.05)
# 
# Aldex_ITSF1_genus.sig<-which(Aldex_ITSF1_genus.all$wi.eBH <=0.05)
# Aldex_ITSF2_genus.sig<-which(Aldex_ITSF2_genus.all$wi.eBH <=0.05)
# Aldex_ITSF3_genus.sig<-which(Aldex_ITSF3_genus.all$wi.eBH <=0.05)
# 
# #Plots of significant families relative abundances
# #Extract significant ASV data from ALDEx2 output
# # 16s
# Aldex_ITSF2_genus.sig.row<-rownames(Aldex_ITSF2_genus.all)[which(Aldex_ITSF2_genus.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_ITSF2_genus.sig.table<-subset(Aldex_ITSF2_genus.all, rownames(Aldex_ITSF2_genus.all) %in% Aldex_ITSF2_genus.sig.row) #Subset significant families
# Aldex_ITSF2_genus.sig.taxon<-subset(taxon.ITS_clean, rownames(taxon.ITS_clean) %in% Aldex_ITSF2_genus.sig.row) #Subset taxon table to get the taxonomy of significant families
# Aldex_ITSF2_genus.sig.table.all<-bind_cols(Aldex_ITSF2_genus.sig.taxon, Aldex_ITSF2_genus.sig.table) #combine tables
# Aldex_ITSF2_genus.sig.table.all$OTU<-rownames(Aldex_ITSF2_genus.sig.table.all)
# 
# Aldex_ITSF3_genus.sig.row<-rownames(Aldex_ITSF3_genus.all)[which(Aldex_ITSF3_genus.all$wi.eBH <=0.05)] #select the rownames of significant families
# Aldex_ITSF3_genus.sig.table<-subset(Aldex_ITSF3_genus.all, rownames(Aldex_ITSF3_genus.all) %in% Aldex_ITSF3_genus.sig.row) #Subset significant families
# Aldex_ITSF3_genus.sig.taxon<-subset(taxon.ITS_clean, rownames(taxon.ITS_clean) %in% Aldex_ITSF3_genus.sig.row) #Subset taxon table to get the taxonomy of significant families
# Aldex_ITSF3_genus.sig.table.all<-bind_cols(Aldex_ITSF3_genus.sig.taxon, Aldex_ITSF3_genus.sig.table) #combine tables
# Aldex_ITSF3_genus.sig.table.all$OTU<-rownames(Aldex_ITSF3_genus.sig.table.all)
# 
# 
# 
# #Effect plots for those ASV with effect size over 1 and less than -1
# Aldex_ITSF2_genus.eff.table.all<-bind_rows(subset(Aldex_ITSF2_genus.sig.table.all, effect >=1),subset(Aldex_ITSF2_genus.sig.table.all, effect <=-1))
# Aldex_ITSF3_genus.eff.table.all<-bind_rows(subset(Aldex_ITSF3_genus.sig.table.all, effect >=1),subset(Aldex_ITSF3_genus.sig.table.all, effect <=-1))
# 
# 
# #Transpose asv table to have ASV in columns
# asv_ITSF2_genus<-t(asv_ITSF2_genus)
# asv_ITSF3_genus<-t(asv_ITSF3_genus)
# 
# #Step 2: Replace zero values before clr transformation. 
# #Uses Count Zero Multiplicative method to replace zeros and outputs pseudo-counts
# asv.n0_ITSF2_genus<-t(cmultRepl(asv_ITSF2_genus, label=0, method="CZM", output="p-counts")) 
# asv.n0_ITSF3_genus<-t(cmultRepl(asv_ITSF3_genus, label=0, method="CZM", output="p-counts"))
# 
# #Note: Check the output to make sure there are no negative numbers. If samples or ASV are sparse, the CZM method will 
# #add a negative number that interferes with the log normalization. If necessary use function below to convert negative values
# #into positives
# asv_n0_ITSF2_genus<-ifelse(asv.n0_ITSF2 < 0, asv.n0_ITSF2*(-1), asv.n0_ITSF2)
# asv_n0_ITSF3_genus<-ifelse(asv.n0_ITSF3 < 0, asv.n0_ITSF3*(-1), asv.n0_ITSF3)
# 
# 
# ## Boxplots of ASVs with effect size grater than 1
# 
# asv.n0.acomp_ITSF2_genus<-as.data.frame(acomp(t(asv_n0_ITSF2_genus)), total=1)
# asv.n0.acomp_ITSF3_genus<-as.data.frame(acomp(t(asv_n0_ITSF3_genus)), total=1)
# 
# #ASV level plots
# 
# #Make Phyloseq object
# phyloseqITSF2_genus <- phyloseq(otu_table(asv.n0.acomp_ITSF2_genus, taxa_are_rows = TRUE), tax_table(taxon_ITS), sample_data(metadata.ITSF2))
# phyloseqITSF3_genus <- phyloseq(otu_table(asv.n0.acomp_ITSF3_genus, taxa_are_rows = TRUE), tax_table(taxon_ITS), sample_data(metadata.ITSF3))
# 
# 
# #ASV level Plots
# #Make long format table from Phyloseq object
# asv_ITSF2_genus_long <- phyloseqITSF2_genus %>%  
#   transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
#   psmelt() %>%  #Melts to long format
#   arrange(desc(Abundance))
# 
# asv_ITSF3_genus_long <- phyloseqITSF3_genus %>%  
#   transform_sample_counts(function(x) {x * 100} ) %>% #Recalculates the relative abundance 
#   psmelt() %>%  #Melts to long format
#   arrange(desc(Abundance))
# 
# 
# #Make vector with ASV names for differential abundant ASVs with effect size over 1 for each pair of comparisons
# asv_ITSF2_genus_effect<-Aldex_ITSF2_genus.eff.table.all$OTU
# asv_ITSF3_genus_effect<-Aldex_ITSF3_genus.eff.table.all$OTU
# 
# 
# #Filter table to obtain only ASVs that with effect size greater than absolute 1
# asv_ITSF2_genus_filter_effect <- filter(asv_ITSF2_genus_long, OTU %in% asv_ITSF2_genus_effect)
# asv_ITSF3_genus_filter_effect <- filter(asv_ITSF3_genus_long, OTU %in% asv_ITSF3_genus_effect)
# 
# 
# #calculate statistics by ASV and Time
# stat_ITSF2_genus<-asv_ITSF2_genus_filter_effect %>%
#   group_by(Time, Genus)%>%
#   summarise(Mean=mean(Abundance, na.rm=TRUE))
# 
# stat_ITSF3_genus<-asv_ITSF3_genus_filter_effect %>%
#   group_by(Time, Genus)%>%
#   summarise(Mean=mean(Abundance, na.rm=TRUE))
# 
# 
# #Reshape table to wide format
# stat_ITSF2_genus_mean<-dcast(stat_ITSF2_genus, Genus~ Time , value.var='Mean')
# stat_ITSF3_genus_mean<-dcast(stat_ITSF3_genus, Genus~ Time , value.var='Mean')
# 
# #Calculate log fold chage - log2 of the ratio of mean RA of first facility by second facility in comparison
# stat_ITSF2_genus_mean$logchange<-log2(stat_ITSF2_genus_mean$Before/stat_ITSF2_genus_mean$After)
# stat_ITSF3_genus_mean$logchange<-log2(stat_ITSF3_genus_mean$Before/stat_ITSF3_genus_mean$After)
# 
# #Add facility identifier - Year with higher RA 
# stat_ITSF2_genus_mean$Time<-ifelse(stat_ITSF2_genus_mean$logchange <0 , "After", "Before")
# stat_ITSF3_genus_mean$Time<-ifelse(stat_ITSF3_genus_mean$logchange <0 , "After", "Before")
# 
# 
# #Plot 
# Logfold_ITSF2_genus<-ggplot(stat_ITSF2_genus_mean, aes(x=logchange, y=reorder(Genus,logchange), fill=Time))+
#   geom_bar(stat='identity', color='black')+
#   geom_text(aes(label=Genus), position = position_stack(vjust=0.5), size=5)+
#   theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=8, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   xlab("Log fold change (log2 Before/After") + 
#   scale_x_continuous(limits=c(-12,12), breaks = c(-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12))+
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="white", color =NA),
#         panel.border = element_rect(color="black", fill=NA)) +
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_manual(values=c("#BB3754","#D6879880"))+
#   #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   ggtitle("Fungi F2 at genus level - Before v After", subtitle = "Differential abundant genera")
# 
# Logfold_ITSF3_genus<-ggplot(stat_ITSF3_genus_mean, aes(x=logchange, y=reorder(Genus,logchange), fill=Time))+
#   geom_bar(stat='identity', color='black')+
#   geom_text(aes(label=Genus), position = position_stack(vjust=0.5), size=5)+
#   theme(axis.title.y = element_blank(), axis.ticks=element_line(color='black'),
#         axis.text=element_text(size=8, color='black')) + 
#   guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
#   xlab("Log fold change (log2 Before/After)") + 
#   scale_x_continuous(limits=c(-12,12), breaks = c(-12,-10,-8,-6,-4,-2,0,2,4,6,8,10,12))+
#   theme(panel.background = element_rect(fill=NA, color =NA),
#         plot.background = element_rect(fill="white", color =NA),
#         panel.border = element_rect(color="black", fill=NA)) +
#   theme(plot.title = element_text(hjust = 0.5, face = 'bold', size=12))+
#   scale_fill_manual(values=c("#FCA50A","#FDC96C80"))+
#   #theme(plot.margin=margin(t=0.2, b=0.2, l=0.2, r=0.2, unit = 'in'))+
#   ggtitle("Fungi F3 genus - Before v After", subtitle = "Differential abundant genera")
# 
# 
# 
# #Combine plots
# Logfold_byTime_genus <- plot_grid(Logfold_ITSF2_genus, Logfold_ITSF3_genus,
#                             ncol=2, nrow=1, labels = c("A","B"), label_size = 20, vjust = 2, hjust = -1.5)
# Logfold_byTime_genus
# ggsave("Aldex_LogFold_byTime_Genus.png", plot=Logfold_byTime_genus, device="png", width=8, height=5, units="in", dpi=600)
# ggsave("Aldex_LogFold_byTime_Genus.svg", plot=Logfold_byTime_genus, device="svg", width=8, height=5, units="in", dpi=600)
