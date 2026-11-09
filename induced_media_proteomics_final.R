library(dplyr)
library(seqinr)
library(devtools)
library(tidyverse)
library(seqinr)
library(splitstackshape)
library(reshape2)
library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(Peptides)
library(Biostrings)
library(purrr)
library(plyr)
library(cluster)
library(factoextra)
library(ggdendro)
library(dendextend)
library(gplots)
library(tidyr)
library(Rmisc)
library(apeglm)
library(pheatmap)
library(limma)
library(edgeR)
library(lsmeans)
library(FSA) 
library(rcompanion)
library(tibble)
library(readxl)

#set working directory
setwd("/Users/Christopher/Dropbox/Proteomics/Induced_media")

#reading in all of the necessary files##################################################
#protein ids that were generated from Rutgers Proteomics
proteomics_number_peptides <- read_xlsx("/Users/Christopher/Dropbox/Proteomics/Induced_media/proteomic_spreadsheets/im_proteomics_jgi_best_models_7april25_long.xlsx")

#DeepLoc output file for predicted protein localization
protein_localization <- read_xlsx("/Users/Christopher/Dropbox/Proteomics/Induced_media/deeploc/jgi_peptidase_s8_deeploc_subcellular_localization.xlsx")
trypsin_localization <- read_xlsx("/Users/Christopher/Dropbox/Proteomics/Induced_media/deeploc/trypsin_deeploc_subcellular_localizations.xlsx")

#metabolomics data
metabolomics <- read.csv("/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/metabolomic_normalized_ion_counts.csv")
  #reading in the sample metadata
meta_data <- read.csv("/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/metabolomics_meta_data.csv")
  #needed for the metabolomics experiment growth curve
induced_media <- read.csv("/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/induced_media_raw_data_file.csv", sep = ";")
  #metabolomics pathway analysis
metabolomic_pathway_analysis <- read.csv("/Users/Christopher Johns/Dropbox/Lith Adsorption/Spreadsheets/induced_media/metabolic_pathway_analysis_updated.csv")
  #ids for pathway analysis
metabo_analyst_ids <- read.csv("/Users/Christopher Johns/Dropbox/Lith Adsorption/Spreadsheets/induced_media/metabolomic_ids.csv")

#transcriptomic datasets
  #ehux transcriptome
ehux_diff_exp <- read.csv("/Users/Christopher/Dropbox/Proteomics/Induced_media/Ehux_TPM_values_w_BLAST.csv")
  #pfam results table
transcriptome_pfam_results <- read.table("/Users/Christopher/Dropbox/Proteomics/Induced_media/pfam_outputs/transcriptome_pfam_annotation", col.names = c("qseqid", "accession", "query_name", "accession", "full_seq_evalue", "full_seq_score", "full_seq_bias", "best_first_domain_evalue", "best_first_domain_score", "best_first_domain_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description of target"), skip = 3)[,-c(2,19)]
  #virus pfam results tables
ehv207_pfam_results <- read.table("ehv_pfam_output/ehv_207_pfam_annotation.txt", col.names = c("qseqid", "accession", "query_name", "accession", "full_seq_evalue", "full_seq_score", "full_seq_bias", "best_first_domain_evalue", "best_first_domain_score", "best_first_domain_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description of target"), skip = 3)[,-c(2,19)]
ehvPS401_pfam_results <- read.csv("ehv_pfam_output/ehv_PS401_pfam_annotation.csv")[,-c(2,19)]
ehv86_pfam_results <- read.csv("ehv_pfam_output/ehv_86_pfam_annotation.csv")[,-c(2,19,20)]
viral_serine_protease_sum <- read.csv("ehv_pfam_output/serine_protease_sum.csv")
 #there were no empty rows without blast hits
empty_blast_hits <- subset(ehux_diff_exp, !complete.cases(ehux_diff_exp$blast_hit))
  #removing columns for the 379 transcriptomic data since they weren't used
ehux_diff_exp_no379 <- subset(ehux_diff_exp, select = -c(X379_C2,X379_C24,X379_V2,X379_V24,log_fold_change,X,X.1,X.2))
  #summing the tpms across all columns
ehux_diff_exp_no379$tpm_sum <- rowSums(ehux_diff_exp_no379[,2:5])
  #removing any zeros that had no mapped reads in the 374 dataset
ehux_diff_exp_no379 <- ehux_diff_exp_no379 %>% filter(tpm_sum >0)
  #there were 85 individual blast hits that mapped to virus genes, 3 were transposons in host genomes and were not counted
virus_genes <- filter(ehux_diff_exp_no379, grepl('virus', blast_hit))

#induced media proteomics#############################################################################################################################################################
#looking at the breakdown of what proteomes the peptides are mapping to
no_ammonium_pie <- proteomics_number_peptides %>% filter(sample == "1") %>% ggplot(aes(x="", y=num_peptides, fill=proteome))+
  geom_bar(width = 1, stat = "identity")+  
  scale_fill_manual(values = c("Emihu" = "#56B4E9", "ehv207" = "#999999"))+
  theme_void()+  
  coord_polar("y", start=0) 


forty_ammonium_pie <- proteomics_number_peptides %>% filter(sample == "3") %>% ggplot(aes(x="", y=num_peptides, fill=proteome))+
  geom_bar(width = 1, stat = "identity")+  
  scale_fill_manual(values = c("Emihu" = "#56B4E9", "ehv207" = "#999999"))+
  theme_void()+  
  coord_polar("y", start=0)

sixty_ammonium_pie <- proteomics_number_peptides %>% filter(sample == "4") %>% ggplot(aes(x="", y=num_peptides, fill=proteome))+
  geom_bar(width = 1, stat = "identity")+
  scale_fill_manual(values = c("Emihu" = "#56B4E9", "ehv207" = "#999999"))+
  theme_void()+  
  coord_polar("y", start=0)

eighty_ammonium_pie <- proteomics_number_peptides %>% filter(sample == "5") %>% ggplot(aes(x="", y=num_peptides, fill=proteome))+
  geom_bar(width = 1, stat = "identity")+ 
  scale_fill_manual(values = c("Emihu" = "#56B4E9", "ehv207" = "#999999"))+
  theme_void()+  
  coord_polar("y", start=0)

#looking at the breakdown of proteins across the different ammonium sulfate fractions, the only fraction to contain viral-derived proteins was the 0% ammonium
ggarrange(no_ammonium_pie, forty_ammonium_pie, sixty_ammonium_pie, eighty_ammonium_pie, nrow = 2, ncol = 2, legend = "bottom", common.legend = TRUE)

#pfam results
jgi_best_genemodel_pfam_results <- read.table("/Users/Christopher/Dropbox/Proteomics/Induced_media/pfam_outputs/pfam_output_proteomics_7april25.txt", col.names = c("qseqid", "accession", "query_name", "accession", "full_seq_evalue", "full_seq_score", "full_seq_bias", "best_first_domain_evalue", "best_first_domain_score", "best_first_domain_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description of target"), skip = 3)[,-c(2,19)]

#taking the lowest evalue
jgi_best_genemodel_pfam_results <- jgi_best_genemodel_pfam_results[jgi_best_genemodel_pfam_results$full_seq_evalue == ave(jgi_best_genemodel_pfam_results$full_seq_evalue, jgi_best_genemodel_pfam_results$qseqid, FUN=min),]

#merge the proteomics results with the pfam output
#renames the qseqid to proteinid
jgi_best_genemodel_pfam_results$protein_id <- jgi_best_genemodel_pfam_results$qseqid
#merges the two files
annotated_proteomics_file <- left_join(proteomics_number_peptides, jgi_best_genemodel_pfam_results)
#replace the NA's with unknown
annotated_proteomics_file <- annotated_proteomics_file %>% replace_na(list(query_name = 'Unknown'))

#plotting the peptide data 
annotated_proteomics_file %>% #removes the viral-associated proteins
  filter(proteome != "ehv207") %>% 
  filter(num_peptides > 1) %>% #removing proteins that have less than two peptides mapped to the protein sequenece
  arrange(desc(num_peptides)) %>%
  mutate(query_name = factor(query_name, levels = unique(query_name))) %>% 
  ggplot(aes(fill=sample_id, y=num_peptides, x=reorder(query_name, -num_peptides, sum)))+ #the minus sign changes it to descending order
  geom_bar(position="stack", stat="identity")+  
  theme_classic2()+
  theme(legend.position = c(0.8,0.7), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 0.9))+  
  scale_y_continuous("Number of unique peptides")+  
  scale_x_discrete("Pfam Query Name")+
  scale_fill_manual(values = c("0_ammonium" = "#0072B2", "40_ammonium" = "#009E73", "60_ammonium" = "#D55E00", "80_ammonium" = "#CC79A7"))

#plotting the protein localization threshold that were found using deeploc 2.0#####################################################################################
#this is for the substilisin
protein_localization$probability <- as.numeric(protein_localization$probability)
protein_localization %>% 
  arrange(protein_id, desc(probability*100)) %>%
  mutate(localizations = factor(localizations, levels = unique(localizations))) %>%
  ggplot(aes(fill=protein_id_short, y=(probability*100), x=localizations))+ #the minus sign changes it to descending order
  geom_bar(stat="identity", position = "dodge")+  
  theme_classic2()+
  theme(legend.position = c(0.8,0.7), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 0.9))+ 
  scale_y_continuous("Localization probability percentages")+
  scale_x_discrete("Subcellular localization")+
  scale_fill_manual(values = c("Emihu_48500016" = "#882255", "Tomato_PR-P69" = "#44AA99", "threshold" = "#D55E00"))

#this is for the trypsin
#pie chart
pie_data <- data.frame(
  protein_id = c("Emihu_1740020","Emihu_840140"),
  number_of_unique_peptides=c(7,4))

colors <- c("#88CCEE","#CC6677")

# Basic piechart
pie_chart <- ggplot(pie_data, aes(x="", y=number_of_unique_peptides, fill=protein_id)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void()+ # remove background, grid, numeric labels
  scale_fill_manual(values=colors)

trypsin_localization$probability <- as.numeric(trypsin_localization$probability)

localizations <- trypsin_localization %>% 
  arrange(protein_id_short, desc(probability*100)) %>%
  mutate(localizations = factor(localizations, levels = unique(localizations))) %>%
  ggplot(aes(fill=protein_id_short, y=(probability*100), x=localizations))+ #the minus sign changes it to descending order
  geom_bar(stat="identity", position = "dodge")+  
  theme_classic2()+
  theme(legend.position = c(0.8,0.7), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 0.9))+ 
  scale_y_continuous("Localization probability percentages", limits = c(0,80))+
  scale_x_discrete("Subcellular localization")+
  scale_fill_manual(values = c("Emihu_1740020" = "#88CCEE", "Emihu_840140" = "#CC6677", "threshold" = "#D55E00"))


ggarrange(pie_chart, localizations, nrow = 2, ncol = 1)

#transcriptome#####################################################################################################################################################################
#read in pfam results of translated protein sequences
#taking the lowest evalue
# pfam_results <- pfam_results[pfam_results$full_seq_evalue == ave(pfam_results$full_seq_evalue, pfam_results$qseqid, FUN=min),]
#there was an _1_1 added to all the Gene ids when I did the protein translation and I need to remove them so they will merge
transcriptome_pfam_results$qseqid <- gsub("_1_1", "", transcriptome_pfam_results$qseqid)

#pfam results from the virus proteomes
ehv207_pfam_results$virus_name <- rep(c("ehv207"))
ehv207_pfam_results$qseqid <- as.factor(ehv207_pfam_results$qseqid)

ehvPS401_pfam_results$virus_name <- rep(c("ehvPS401"))
ehvPS401_pfam_results$qseqid <- as.factor(ehvPS401_pfam_results$qseqid)

ehv86_pfam_results$virus_name <- rep(c("ehv86"))
ehv86_pfam_results$qseqid <- as.factor(ehv86_pfam_results$qseqid)

#bind datasets
viral_proteomes <- bind_rows(ehv207_pfam_results,ehv86_pfam_results,ehvPS401_pfam_results)
viral_proteomes_serine_proteaes <- subset(viral_proteomes, str_detect(query_name, 'Trypsin') | str_detect(query_name, 'Peptidase_S') | str_detect(query_name, 'CLP_protease'))

#plot for figure
viral_serine_protease_sum %>%
  filter(virus != "ehvps401") %>%
  ggplot(aes(fill=protease, x=virus, y=number))+ 
  geom_bar(stat="identity", color ="black")+
  theme_classic2()+
  theme(legend.position = c(0.8,0.7), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  scale_y_continuous("Number of Serine Proteases")+  
  scale_x_discrete("Virus", labels = c("EhV-207","EhV-86","EhV-PS401"))+  
  scale_fill_manual(values = c("CLP Protease" = "#004949", "Trypsin" = "#006ddb", "Peptidase S49" = "#920000"))

#to look at normalized values across samples to account for any differences in the total number of normalized counts, tpm, between samples
ehux_diff_exp$psuedo_reference_sample_24h <- sqrt((ehux_diff_exp$X374_V24+0.00001)*(ehux_diff_exp$X374_C24+0.00001))
ehux_diff_exp$psuedo_reference_sample_2h <- sqrt((ehux_diff_exp$X374_V2+0.00001)*(ehux_diff_exp$X374_C2+0.00001))
#calculating the ratio of each sample to the reference
#24h
ehux_diff_exp$ratio_24h_control <- (ehux_diff_exp$X374_C24/ehux_diff_exp$psuedo_reference_sample_24h)
ehux_diff_exp$ratio_24h_infected <- (ehux_diff_exp$X374_V24/ehux_diff_exp$psuedo_reference_sample_24h)
#2h
ehux_diff_exp$ratio_2h_control <- (ehux_diff_exp$X374_C2/ehux_diff_exp$psuedo_reference_sample_2h)
ehux_diff_exp$ratio_2h_infected <- (ehux_diff_exp$X374_V2/ehux_diff_exp$psuedo_reference_sample_2h)
#calculating the normalization factor for each sample
normalization_factor_control24 <- median(ehux_diff_exp$ratio_24h_control) 
normalization_factor_infected24 <- median(ehux_diff_exp$ratio_24h_infected) 
normalization_factor_control2 <- median(ehux_diff_exp$ratio_2h_control) 
normalization_factor_infected2 <- median(ehux_diff_exp$ratio_2h_infected)
#size factors are pretty similar so I'll go ahead and just use the TPM values for fold change

#adding small number to the TPM
ehux_diff_exp$fold_change_24 <- (ehux_diff_exp$X374_V24+0.00001)/(ehux_diff_exp$X374_C24+0.00001)
ehux_diff_exp$fold_change_2 <- (ehux_diff_exp$X374_V2+0.00001)/(ehux_diff_exp$X374_C2+0.00001)
ehux_diff_exp$fold_change_2_379 <- (ehux_diff_exp$X379_V2+0.00001)/(ehux_diff_exp$X379_C2+0.00001)
ehux_diff_exp$fold_change_24_379 <- (ehux_diff_exp$X379_V24+0.00001)/(ehux_diff_exp$X379_C24+0.00001)

#calculation log2(tpm + 0.00001)
ehux_diff_exp$log2_fold_change_2 <- log2(ehux_diff_exp$fold_change_2)
ehux_diff_exp$log2_fold_change_24 <- log2(ehux_diff_exp$fold_change_24)
ehux_diff_exp$log2_fold_change_2_379 <- log2(ehux_diff_exp$fold_change_2_379)
ehux_diff_exp$log2_fold_change_24_379 <- log2(ehux_diff_exp$fold_change_24_379)

only374 <- filter(ehux_diff_exp, X374_C2 == "0" & X374_C24 == "0" & X374_V2 == "0" & X374_V24 == "0")

#dropping unecessary columns
ehux_diff_exp <- within(ehux_diff_exp, rm(log_fold_change, X, X.1, X.2))
ehux_diff_exp$Gene <- gsub("_1", "", ehux_diff_exp$Gene)
ehux_diff_exp$Gene <- as.character(ehux_diff_exp$Gene)

transcriptome_pfam_results$Gene<- transcriptome_pfam_results$qseqid
transcriptome_pfam_results$Gene <- as.character(transcriptome_pfam_results$Gene)

ehux_diff_exp <- left_join(ehux_diff_exp, transcriptome_pfam_results, by = "Gene")

#pulling out all of the serine proteases
ehux_diff_exp_serine_proteases <- subset(ehux_diff_exp, str_detect(query_name, 'Peptidase_S') | str_detect(query_name, 'Trypsin') | str_detect(query_name, 'CLP_protease'))
ehux_diff_exp_serine_proteases <- filter(ehux_diff_exp_serine_proteases, log2_fold_change_24 != "0")

ehux_diff_exp_subtilisin_proteases <- subset(ehux_diff_exp_serine_proteases, str_detect(query_name, 'Peptidase_S8'))

ehux_diff_exp_serine_proteases <- ehux_diff_exp_serine_proteases[ehux_diff_exp_serine_proteases$full_seq_evalue == ave(ehux_diff_exp_serine_proteases$full_seq_evalue, ehux_diff_exp_serine_proteases$qseqid, FUN=min),]

ehux_diff_exp_serine_proteases <- ehux_diff_exp_serine_proteases[ehux_diff_exp_serine_proteases$best_first_domain_evalue == ave(ehux_diff_exp_serine_proteases$best_first_domain_evalue, ehux_diff_exp_serine_proteases$qseqid, FUN=min),]

ehux_diff_exp_serine_proteases$diff_exp <- ifelse(ehux_diff_exp_serine_proteases$log2_fold_change_2<0, "down",
                                                  ifelse(ehux_diff_exp_serine_proteases$log2_fold_change_2>0, "up",
                                                         ifelse(ehux_diff_exp_serine_proteases$log2_fold_change_2 == 0, "no expression", "")))

ehux_diff_exp_serine_proteases$diff_exp_24 <- ifelse(ehux_diff_exp_serine_proteases$log2_fold_change_24<0, "down",
                                                  ifelse(ehux_diff_exp_serine_proteases$log2_fold_change_24>0, "up",
                                                         ifelse(ehux_diff_exp_serine_proteases$log2_fold_change_24 == 0, "no expression", "")))
#cleaning up the protease names
ehux_diff_exp_serine_proteases$step1 <- gsub("_2", "", ehux_diff_exp_serine_proteases$query_name)
ehux_diff_exp_serine_proteases$protease_name <- gsub("_N", "", ehux_diff_exp_serine_proteases$step1)

#saving the serine protease protein sequences for the NAAMES MetaT
# ehux_metaT_prot_seq <- left_join(ehux_diff_exp_serine_proteases, trans_prot_seq, by = "Gene")
# ehux_metaT_prot_seq <- subset(ehux_metaT_prot_seq, select = c("Gene","prot_seq"))

#plot for fig
serine_protease_24h <- ehux_diff_exp_serine_proteases %>%
  mutate(Gene = factor(Gene, levels = unique(Gene))) %>% 
  ggplot(aes(fill=diff_exp, y=log2_fold_change_24, x=reorder(Gene, -log2_fold_change_24)))+ #the minus sign changes it to descending order
  geom_abline(slope=0, intercept=0,  col = "black")+
  geom_abline(slope=0, intercept=5,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=-5,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=10,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=-10,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=15,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=-15,  col = "grey",lty=2)+
  geom_bar(stat="identity", color ="black")+
  theme_classic2()+
  theme(legend.position = c(0.8,0.7), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),  
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_continuous("Fold change; log2")+  
  scale_x_discrete("Gene")+  
  scale_fill_manual(values = c("down" = "#00AFBB", "up" = "#FC4E07", "no expression" = "#E7B800"))

#source data
# write.csv(ehux_diff_exp_serine_proteases, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig1b.csv")

#splitting the data set into host and viral derived proteins for clustering
virus_diff_exp_serine_protease <- subset(ehux_diff_exp_serine_proteases, str_detect(blast_hit, 'virus'))
ehux_diff_exp_serine_proteases <- ehux_diff_exp_serine_proteases %>% filter(blast_hit != ">CAZ69462 putative ATP-dependent protease proteolytic subunit [Emiliania huxleyi virus 99B1].")

#looking at the distribution of fold changes
hist(ehux_diff_exp_serine_proteases$log2_fold_change_2, breaks = 3)
shapiro.test(ehux_diff_exp_serine_proteases$log2_fold_change_2)
shapiro.test(ehux_diff_exp_serine_proteases$log2_fold_change_24)

#clustering data
clustering <- subset(ehux_diff_exp_serine_proteases, select = c("Gene", "protease_name", "log2_fold_change_2", "log2_fold_change_24"))
clustering$ID <- as.factor(paste(clustering$protease_name, clustering$Gene, sep="-"))

#need rownames for the cluster analysis
rownames(clustering) = clustering$ID
#subset the fold changes
clustering <- subset(clustering, select = c("log2_fold_change_2","log2_fold_change_24"))
#estimating the number of clusters using the elbow method
fviz_nbclust(clustering, kmeans, method = "wss", nboot = 500, k.max = 15)
#hierachal clustering for the dendrogram
d <- dist(clustering, method = "euclidean") 
#perform hierarchical clustering using Ward's method
final_clust <- hclust(d, method = "ward.D2" )
#cut the dendrogram into 6 clusters based on the elbow plot
groups <- cutree(final_clust, k=6)
plot(final_clust)

#looking at the dendrogram
dendro <- as.dendrogram(final_clust)
dendro_plot <- ggdendrogram(data = dendro, rotate = TRUE)
dendro_plot <- dendro_plot + theme(axis.text.y = element_text(size = 6))

#find number of observations in each cluster
table(groups)

#add the group to the data frame
final_data <- cbind(ehux_diff_exp_serine_proteases, cluster = groups)

#plotting as a heatmap
heat_map_data <- subset(final_data, select = c("Gene","blast_hit","query_name","log2_fold_change_2","log2_fold_change_24","cluster"))
heat_map_data_virus <- subset(virus_diff_exp_serine_protease, select = c("Gene","blast_hit","query_name","log2_fold_change_2","log2_fold_change_24"))
heat_map_data <- reshape2::melt(heat_map_data, id.vars = c("Gene","blast_hit","query_name","cluster"), value.name = "log2_fold_change")
heat_map_data_virus <- reshape2::melt(heat_map_data_virus, id.vars = c("Gene","blast_hit","query_name"), value.name = "log2_fold_change")

heat_map_data$step1 <- gsub("_2", "", heat_map_data$query_name)
heat_map_data_virus$step1 <- gsub("_2", "", heat_map_data_virus$query_name)

heat_map_data$protease_name <- gsub("_N", "", heat_map_data$step1)
heat_map_data_virus$protease_name <- gsub("_N", "", heat_map_data_virus$step1)

heat_map_data_rescale <- heat_map_data

#binning the values for visualization
heat_map_data_rescale$range <- cut(heat_map_data_rescale$log2_fold_change, breaks = c(-20,-15,-10,-5,-1,-0.5,-0.25,0,0.25,0.5,1,5,10,15,20,25,30))
heat_map_data_virus$range <- cut(heat_map_data_virus$log2_fold_change, breaks = c(-20,-15,-10,-5,-1,-0.5,-0.25,0,0.25,0.5,1,5,10,15,20,25,30))

#creating new ids
heat_map_data_rescale$ID <- as.factor(paste(heat_map_data$protease_name, heat_map_data$Gene, sep="-"))
heat_map_data_virus$ID <- as.factor(paste(heat_map_data_virus$protease_name, heat_map_data_virus$Gene, sep="-"))

#join the datasets together
heat_map_data_rescale <- left_join(heat_map_data_rescale, ehux_diff_exp_serine_proteases, by = "Gene")

#generating the final plots for the manuscript
#order the rows by the dendrogram
dendro_order <- order.dendrogram(dendro)
heat_map_data_rescale$ID <- factor(x = heat_map_data_rescale$ID,
                               levels = heat_map_data_rescale$ID[dendro_order], 
                               ordered = TRUE)

#e hux-derived heatmap
heatmap_plot <- heat_map_data_rescale %>% ggplot(aes(x = variable, y = ID))+
  geom_tile(aes(fill = range, colour = "black"))  + labs (fill= expression(log[2]~"Fold change")) +
  scale_fill_manual(breaks = levels(heat_map_data_rescale$range),
                    values = c("#33608C", "#626496", "#8867A1", "#AC6A9F", "#CD6E94","#E87687","#F09574","#F5B35E","#F4A652","#F08B4B","#EA7246","#E06046","#D64E45","#CA3B43","#BE2641","#B81840"))+
  #geom_text(size=4, color="white", aes(label=formatC(E_all, format = "e", digits = 1))) +
  # theme_Publication2() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none", 
        legend.direction = "horizontal", legend.key.height=unit(1.5,"line"), legend.key.width = unit (2.5, "line"), plot.background=element_blank())

#this was pulled out since it is EhV-derived
heat_map_data_virus %>% ggplot(aes(x = variable, y = ID))+
  geom_tile(aes(fill = range, colour = "black"))  + labs (fill= expression(log[2]~"Fold change")) +
  scale_fill_manual(breaks = levels(heat_map_data_rescale$range),
                    values = c("#33608C", "#626496", "#8867A1", "#AC6A9F", "#CD6E94","#E87687","#F09574","#F5B35E","#F4A652","#F08B4B","#EA7246","#E06046","#D64E45","#CA3B43","#BE2641","#B81840"))+
  #geom_text(size=4, color="white", aes(label=formatC(E_all, format = "e", digits = 1))) +
  # theme_Publication2() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        legend.position = "none", 
        legend.direction = "horizontal", legend.key.height=unit(1.5,"line"), legend.key.width = unit (2.5, "line"), plot.background=element_blank())

grid.newpage()
print(heatmap_plot, 
      vp = viewport(x = 0.3, y = 0.5, width = 0.5, height = 1.0))
print(dendro_plot, 
      vp = viewport(x = 0.70, y = 0.49, width = 0.3, height = 1.12))

#just the subtilisin figure
ehux_diff_exp_subtilisin <- filter(ehux_diff_exp_serine_proteases, str_detect(query_name, 'Peptidase_S'))
ehux_diff_exp_subtilisin <- filter(ehux_diff_exp_serine_proteases, str_detect(query_name, 'Peptidase_S8'))

#plotting data
subtilisin_24h <- ehux_diff_exp_subtilisin %>%
  mutate(Gene = factor(Gene, levels = unique(Gene))) %>% 
  ggplot(aes(fill=diff_exp, y=log2_fold_change_24, x=reorder(Gene, -log2_fold_change_24)))+ #the minus sign changes it to descending order
  geom_abline(slope=0, intercept=0,  col = "black")+
  geom_abline(slope=0, intercept=0.5,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=-0.5,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=1,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=-1,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=1.5,  col = "grey",lty=2)+
  geom_abline(slope=0, intercept=-1.5,  col = "grey",lty=2)+
  geom_bar(stat="identity", color ="black")+  
  theme_classic2()+
  theme(legend.position = c(0.8,0.7), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_y_continuous("Fold change; log2")+  
  scale_x_discrete("Gene ID")+
  scale_fill_manual(values = c("down" = "#00AFBB", "up" = "#FC4E07", "no expression" = "#E7B800"))

#source data
# write.csv(ehux_diff_exp_subtilisin, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig1e.csv")

#analyzing the metabolomic data##########################################################################################################
#reading in the cell growth data from the metabolomics experiment
metabolomics_growth <- subset(induced_media, experiment.name == "metabolomics_exp")

#summarizing the growth data
metabolomics_growth_sum <- metabolomics_growth %>% summarySE(measurevar = "cell.concentration",  
                                     groupvars = c("hour","treatment"))

#growth curve for the metabolimcs experiment
growth_curve <- metabolomics_growth_sum %>%  
  # filter(hour < 144) %>%
ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 12, angle = 45),
        axis.text.y = element_text(size = 12),  
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)",  
                     breaks=c(0,24,48,72,96,120,144,168,192))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment", 
                      values = c("control" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "evs" = "#D55E00"), guide = "none")

#processing the metabolites
#adding up all of the normalized ion counts in each row
metabolomics$ion_sum <- rowSums(metabolomics[ , c(2:19)])

#removing the standards from the data
metabolomics <- metabolomics[!grepl("PosIS",metabolomics$compound),]
metabolomics <- metabolomics[!grepl("NegIS",metabolomics$compound),]

#selecting the polarity with the highest normalized ion counts
metabolomics_dup_rem <- metabolomics[metabolomics$ion_sum == ave(metabolomics$ion_sum, metabolomics$compound, FUN=max),]
metabolomics_dup_rem <- metabolomics_dup_rem %>% filter(compound != "Glycerophosphocholine (+HAc)")
# write.csv(metabolomics_dup_rem, file = "/Users/Christopher Johns/Dropbox/Lith Adsorption/Spreadsheets/induced_media/metabolome.csv")

#for PCA- this removes the columns with the ion sum and the polarity
metabolomics_pca <- subset(metabolomics_dup_rem, select = c(1:19))
#removes the column with metabolite names
CleanInput <- metabolomics_pca[,c(2,3:(ncol(metabolomics_pca)))]
#coverts to a data matrix
DataMatrix <- data.matrix(CleanInput[,1:ncol(CleanInput)])
#adds the rownames back in
row.names(DataMatrix) <- metabolomics_pca[,1]
#pca of the scaled transposed ion counts using prcomp
res.pca <- prcomp(t(DataMatrix), scale = TRUE)
#visualizing the percentage of explained variance 
fviz_eig(res.pca)
#initial look at the difference samples with two dimensions
fviz_pca_ind(res.pca, col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping
#prepping dataset for plotting
PCi<-data.frame(res.pca$x,Treatment=c(rep("Cells-Ctrl",3),rep("Cells-+Protein",3),rep("Cells-+Vesicles",3)))
#transferring row names to columns
PCi <- rownames_to_column(PCi, "sample")
#writing the PCA table to a file
write.csv(PCi,file="PCi Table.csv")
#making the final plot
pca_plot <- PCi %>% dplyr::slice(1:9) %>%
  ggplot(aes(x=PC1,y=PC2,col=Treatment,label=sample)) +
  geom_point(size=5)+ 
  scale_x_continuous(name = "PC1 (58.9%)") +
  scale_y_continuous(name = "PC2 (22%)") +
  theme_classic2()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),  
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  scale_color_manual(values =c("Cells-Ctrl" = "#FFCC00",
                               "Cells-+Protein" = "#FF3300",
                               "Cells-+Vesicles" = "#D55E00"), guide = "none")

ggarrange(growth_curve, pca_plot, nrow = 2, ncol = 2, align = 'h')

#separating into cell and supernatant
metabolomics_cells <- subset(metabolomics_dup_rem, select = c(1:10,20))
# metabolomics_sup <- subset(metabolomics_dup_rem, select = c(1,11:20))

row.names(metabolomics_cells) <- metabolomics_cells[,1]
metabolomics_cells_row_names <- subset(metabolomics_cells,  select = c(2:10))

#using limma for the differential accumulation of metabolites
#I first normalized the ion counts using quartile  normalization, 
#It assumes that the global differences in distributions between samples 
#are all due to technical variation. Any remaining differences are likely actual biological effects.
#Also normalizes the means across samples
#this was followed by a log2 transformation prior to differential expression analysis using limma

#protein vs control
Group <- rep(c("control","protein"),each=3)
design <- model.matrix(~Group)

metabolomic_protein <- subset(metabolomics_cells_row_names, select = c(1:6))

y <- normalizeQuantiles(metabolomic_protein)
y = log2(y+0.001)
boxplot(y)

protein_fit <- lmFit(y, design)
protein_fit <- eBayes(protein_fit)
protein_results <- topTable(protein_fit,coef=2,number=nrow(protein_fit))
head(protein_results)

#vesicles vs control
Group <- rep(c("control","vesicles"),each=3)
design <- model.matrix(~Group)

metabolomic_vesicles <- subset(metabolomics_cells_row_names, select = c(1:3,7:9))

y <- normalizeQuantiles(metabolomic_vesicles)
y = log2(y+0.001) 
boxplot(y)

vesicle_fit <- lmFit(y, design)
vesicle_fit <- eBayes(vesicle_fit)
summary(decideTests(vesicle_fit))
vesicle_results <- topTable(vesicle_fit,coef=2,number=nrow(vesicle_fit))
head(vesicle_results)

#protein vs vesicles
Group <- rep(c("protein","vesicles"),each=3)
design <- model.matrix(~Group)

metabolomic_vesicles_protein <- subset(metabolomics_cells_row_names, select = c(4:9))

y <- normalizeQuantiles(metabolomic_vesicles_protein)
y = log2(y+0.001)  
boxplot(y)

vesicle_protein_fit <- lmFit(y, design)
vesicle_protein_fit <- eBayes(vesicle_protein_fit)
vesicle_protein_results <- topTable(vesicle_protein_fit,coef=2,number=nrow(vesicle_protein_fit))
head(vesicle_protein_results)

#generating the volcano plots
protein_volcano <- protein_results %>% rownames_to_column("compound") %>% mutate(Color = ifelse(P.Value<0.05 & (abs(logFC)) > 1, "yes", "no")) %>%
  ggplot(aes(x = logFC, y = -(log10(P.Value)), colour = Color, label = compound))+
  geom_point(size = 3)+  
  theme_bw()+  
  theme_classic2()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),  
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_hline(yintercept=-(log10(0.05)),linetype=2)+ 
  geom_vline(xintercept=1,linetype=2)+
  geom_vline(xintercept=-1,linetype=2)+
  scale_color_manual(values = c("yes" = "red", "no" = "blue"), guide = "none")+
  geom_text(aes(label = ifelse(P.Value<0.05 & (abs(logFC)) > 1 ,as.character(compound),"")), size = 3, hjust = 0.5, vjust = -0.5)


vesicle_volcano <- vesicle_results %>% rownames_to_column("compound") %>% mutate(Color = ifelse(P.Value<0.05 & (abs(logFC)) > 1 , "yes", "no")) %>%
  ggplot(aes(x = logFC, y = -(log10(P.Value)), colour = Color, label = compound))+
  geom_point(size = 3)+  
  theme_bw()+  
  theme_classic2()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),  
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_hline(yintercept=-(log10(0.05)),linetype=2)+ 
  geom_vline(xintercept=1,linetype=2)+
  geom_vline(xintercept=-1,linetype=2)+
  scale_color_manual(values = c("yes" = "red", "no" = "blue"), guide = "none")+
  geom_text(aes(label = ifelse(P.Value<0.05 & (abs(logFC)) > 1 ,as.character(compound),"")), size = 3, hjust = 0.5, vjust = -0.5)


vesicle_protein_results %>% rownames_to_column("compound") %>% mutate(Color = ifelse(P.Value < 0.05, "yes", "no")) %>%
  ggplot(aes(x = logFC, y = -(log10(P.Value)), colour = Color, label = compound))+
  geom_point(size = 3)+  
  theme_bw()+  
  theme_classic2()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),  
        axis.title=element_text(size=12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_hline(yintercept=-(log10(0.05)),linetype=2)+ 
  geom_vline(xintercept=1,linetype=2)+
  geom_vline(xintercept=-1,linetype=2)+
  scale_color_manual(values = c("yes" = "red", "no" = "blue"), guide = "none")+
  geom_text(aes(label = ifelse(P.Value<0.05 & (abs(logFC)) > 0.5 ,as.character(compound),"")), size = 3, hjust = 0.5, vjust = -0.5)

#plotting everything together
ggarrange(growth_curve, pca_plot, vesicle_volcano, protein_volcano, nrow = 2, ncol = 2, align = 'hv')

#metabolic pathway analysis for overexpressed metabolites
vesicle_pathway_up <- metabolomic_pathway_analysis %>% filter(treatment == "vesicles_vs_control" & exp_profile == "up") %>%
ggplot(aes(x= impact, y = log_p, size = impact, color = log_p)) +
  geom_point() +
  theme_bw()+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  scale_y_continuous()+
  scale_size(range = c(3,15), guide = "none")+
  scale_color_continuous(low = "blue", high = "red")+
  geom_text(aes(label = pathway_name), size = 3, hjust = -0.2, vjust = -0.8)

vesicle_pathway_down <- metabolomic_pathway_analysis %>% filter(treatment == "vesicles_vs_control" & exp_profile == "down") %>%
  ggplot(aes(x= impact, y = log_p, size = impact, color = log_p)) +
  geom_point() +
  theme_bw()+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  scale_y_continuous()+
  scale_size(range = c(3,15), guide = "none")+
  scale_color_continuous(low = "blue", high = "red")+
  geom_text(aes(label = pathway_name), size = 3, hjust = -0.2, vjust = -0.8)

vesicle_pathway <- metabolomic_pathway_analysis %>% filter(treatment == "vesicles_vs_control") %>%
  ggplot(aes(x= impact, y = log_p., size = impact, color = log.p.)) +
  geom_point() +
  theme_bw()+
  theme_classic2()+
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  scale_size(range = c(3,15), guide = "none")+
  scale_color_continuous(low = "blue", high = "red", limits = c(0,3))+
  geom_text(aes(label = pathway_name), size = 3, hjust = -0.2, vjust = -0.8)

ggarrange(vesicle_pathway_up, vesicle_pathway_down, nrow = 2, ncol =2, legend = "bottom" )