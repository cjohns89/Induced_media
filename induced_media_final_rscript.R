library(ggplot2)
library(Rmisc)
library(dplyr)
library(gridExtra)
library(gtable)
library(ggpubr)
library(grid)
library("cowplot")
library(reshape2)
library(data.table)
library(plyr)
library(expss)
library(multcomp)
library(multcompView)
library(lawstat)
library(FSA)
library(rcompanion)
library(broom)
library(lemon)
library(doBy)
library(drda)
library(car)
library(patchwork)
library(ggplotify)
library(agricolae)
library(rstatix)
library(tidyverse)
library(agricolae)
library(userfriendlyscience)
library(ARTool)
library(patchwork)
library(ggplotify)
library(readxl)

#set directory
setwd()

#read in necessary files
induced_media <- read.csv("induced_media_raw_data_file.csv", sep = ";")
tff_protein_concentrations <- read.csv("20212710_tff_im_concentrate.csv")
subtilisin_exp <- read.csv("naames_data/station_4_subtilisin_expression.csv")
subtilisin_exp_updated <- read.csv("naames_data/Ehuxchangedsub1023.csv")
cell_viability_anti_viral <- read.csv("20221410_cell_viability_assay.csv")
protease_activity <- read_excel("protease_activity_measurements.xlsx")
antibiotic_bacterial_concentrations <- read.csv("antibiotic_treatment_bacteria_concentrations.csv")

#the impact of induced media on cell physiology######################################################################################
IM <- induced_media %>% filter(treatment == "control im" | treatment == "viral im") %>% filter(experiment.name != "IM")
IM$hour <- as.numeric(IM$hour)
IM$plotID <- as.factor(paste(IM$treatment, IM$liths, sep="-"))

calcified_sum <- IM %>%  
  filter(hour < 121) %>%
  filter(experiment.name != "haptophyte_screen") %>%
  filter(experiment.name != "haptophyte_screen_2") %>%
  summarySE(measurevar = "calcified", groupvars = c("hour", "treatment", "liths" , "plotID"), na.rm = TRUE)


calcified <- calcified_sum %>%
  arrange(calcified) %>% 
  mutate(plotID = factor(plotID, levels=c("control im-control",  
                                          "control im-cells+liths",
                                          "viral im-control",  
                                          "viral im-cells+liths"))) %>%  
  ggplot(aes(x=hour, y=calcified, colour=plotID, shape = liths))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=calcified-se, ymax=calcified+se), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+  
  scale_y_continuous(name = "Percent calcified cells", limits = c(0,70), breaks = c(0,10,20,30,40,50,60,70))+
  scale_colour_manual("plotID",
                      # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                      values = c("control im-control" = "#006699",
                                 "viral im-control" = "#FF3300",
                                 "control im-cells+liths" = "#006699",  
                                 "viral im-cells+liths" = "#FF3300"))+  
  scale_shape_manual(values = c(1,19,1,19))

#write source data file  
# write.csv(calcified_sum, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs4c.csv")

IM_cells <- IM  %>% filter(experiment.name != "hapto_screen") %>% filter(experiment.name != "hapto_screen_2") %>% 
  summarySE(measurevar = "cell.concentration",  
                           groupvars = c("hour","treatment","liths","plotID"))

#this is for summary purposes to find the number of data points for each treatment
IM %>% filter(experiment.name != "IM") %>% filter(experiment.name != "hapto_screen") %>% filter(experiment.name != "hapto_screen_2") %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("treatment","liths","plotID"))

IM %>% filter(experiment.name != "IM") %>% filter(experiment.name != "hapto_screen") %>% filter(experiment.name != "hapto_screen_2") %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("experiment.name", "treatment","liths","plotID"))

IM_liths <- IM %>% filter(experiment.name != "IM") %>% filter(experiment.name != "hapto_screen") %>% filter(experiment.name != "hapto_screen_2") %>%
  summarySE(measurevar = "lith.concentration",  
                      groupvars = c("hour","treatment","liths","plotID"))
IM_liths %>% filter(liths != "control")

IM_vesicle <- induced_media %>% filter(experiment.name == "proliths_vesicle") %>%
  summarySE(measurevar = "vesicle.concentration",  
            groupvars = c("hour","treatment","liths"))
IM_vesicle$plotID <- as.factor(paste(IM_vesicle$treatment, IM_vesicle$liths, sep="-"))


#plotting the growth curve 
cells <- IM_cells %>% filter(hour < 144) %>%
  arrange(cell.concentration) %>% 
  mutate(plotID = factor(plotID, levels=c("control im-control",  
                                          "control im-cells+liths",
                                          "viral im-control",  
                                          "viral im-cells+liths"))) %>%  
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=plotID, shape = liths))+  
  geom_point(size=5)+theme_bw()+   
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
scale_colour_manual("plotID",
                    # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                    values = c("control im-control" = "#006699",
                               "viral im-control" = "#FF3300",
                               "control im-cells+liths" = "#006699",  
                               "viral im-cells+liths" = "#FF3300"))+  
  scale_shape_manual(values = c(1,19,1,19))

# write.csv(IM_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig2a.csv")

liths <- IM_liths %>% filter(hour < 144) %>%
  arrange(lith.concentration) %>%
  mutate(plotID = factor(plotID, levels=c("control im-cells+liths",  
                                          "viral im-cells+liths"))) %>%  
  ggplot(aes(x=hour, y=log10(lith.concentration), colour=plotID, shape = liths))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(lith.concentration-se), ymax=log10(lith.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120))+  
  scale_y_continuous(name = bquote(atop('Lith concentration; log10' , "coccoliths mL"^-1)))+
  scale_colour_manual("plotID",
                      labels = c("OM + coccoliths", "IM + coccoliths"),
                      values = c("control im-cells+liths" = "#006699",  
                                 "viral im-cells+liths" = "#FF3300"))+  
  scale_shape_manual(values = c(1,1))

#write source data file
# write.csv(IM_liths, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs4d.csv")

ggarrange(cells, calcified, align = 'h', label.x = 0, 
          nrow = 2, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = TRUE, legend = "bottom")

#sytox
IM_sytox <- IM %>% filter(experiment.name == "IM stress") %>% summarySE(measurevar = "normalized.sytox", groupvars = c("hour","treatment","liths","plotID"))

sytox <- IM_sytox %>% filter(hour < 144) %>% 
  arrange(normalized.sytox) %>%
  mutate(plotID = factor(plotID, levels=c("control im-control",  
                                          "control im-cells+liths",
                                          "viral im-control",  
                                          "viral im-cells+liths"))) %>% 
  ggplot(aes(x=hour, y=normalized.sytox, colour=plotID, shape = liths))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=normalized.sytox-se, ymax=normalized.sytox+se), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+
  # geom_hline(yintercept = sytox_mean, color="black", 
  #             linetype="dashed", size=1.5)+
  scale_x_continuous(name = "Hour", limits = c(-10,120), breaks=c(0,24,48,72,96,120))+  
  scale_colour_manual("plotID",
                      # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                      values = c("control im-control" = "#006699",
                                 "viral im-control" = "#FF3300",
                                 "control im-cells+liths" = "#006699",  
                                 "viral im-cells+liths" = "#FF3300"))+  
  scale_shape_manual(values = c(1,19,1,19))+
  scale_y_continuous(name = "Percentage of dead cells", limits = c(0,20))

write.csv(IM_sytox, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs2c.csv")

#ros

IM_ros <- IM %>% filter(experiment.name == "IM stress") %>% summarySE(measurevar = "normalized.ros", groupvars = c("hour","treatment","liths","plotID"))

ros <- IM_ros %>% filter(hour < 144) %>% 
  arrange(normalized.ros) %>%
  mutate(plotID = factor(plotID, levels=c("control im-control",  
                                          "control im-cells+liths",
                                          "viral im-control",  
                                          "viral im-cells+liths"))) %>% 
  ggplot(aes(x=hour, y=normalized.ros, colour=plotID, shape = liths))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+   
  geom_errorbar(aes(ymin=normalized.ros-se, ymax=normalized.ros+se), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+  
  scale_x_continuous(name = "Hour", limits = c(-10,120), breaks=c(0,24,48,72,96,120))+  
  scale_colour_manual("plotID",
                      # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                      values = c("control im-control" = "#006699",
                                 "viral im-control" = "#FF3300",
                                 "control im-cells+liths" = "#006699",  
                                 "viral im-cells+liths" = "#FF3300"))+  
  scale_shape_manual(values = c(1,19,1,19))+
  scale_y_continuous(name = "Reactive oxygen stress \n(median fluorescence units x 1000)", breaks = c(0,1000,2000), labels = c(0,1,2))

# write.csv(IM_ros, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig2d.csv")

#lysotracker
IM$normalized.lyso <- IM$lyso - IM$background_per_lyso_mito
IM$normalized.lyso_flu <- IM$lyso_flu - IM$background_flu_lyso_mito

IM_lyso <- IM %>% filter(experiment.name == "im_lyso_mito") %>% summarySE(measurevar = "normalized.lyso", groupvars = c("hour","treatment","liths","plotID"))

lyso <- IM_lyso %>% filter(hour < 144) %>%
  arrange(normalized.lyso) %>%
  mutate(plotID = factor(plotID, levels=c("control im-control",  
                                          "control im-cells+liths",
                                          "viral im-control",  
                                          "viral im-cells+liths"))) %>% 
  ggplot(aes(x=hour, y=normalized.lyso, colour=plotID, shape = liths))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=normalized.lyso-se, ymax=normalized.lyso+se), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+  
  scale_x_continuous(name = "Hour", limits = c(-10,130), breaks=c(0,24,48,72,96,120))+  
  scale_colour_manual("plotID",
                      values = c("control im-control" = "#006699",
                                 "viral im-control" = "#FF3300",
                                 "control im-cells+liths" = "#006699",  
                                 "viral im-cells+liths" = "#FF3300"))+  
  scale_shape_manual(values = c(1,19,1,19))+
  scale_y_continuous(name = "Percentage of cells\nwith increased lysosomes")

#source data
# write.csv(IM_lyso, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs4a.csv")

#mitotracker
IM_mito <- IM %>% filter(experiment.name == "im_lyso_mito") %>% summarySE(measurevar = "high_mito", groupvars = c("hour","treatment","liths","plotID"))

IM$mito_high <- as.numeric(IM$mito_high)

mito <- IM_mito %>% filter(hour < 144) %>% filter(hour != "72") %>%
  arrange(high_mito) %>%
  mutate(plotID = factor(plotID, levels=c("control im-control",  
                                          "control im-cells+liths",
                                          "viral im-control",  
                                          "viral im-cells+liths"))) %>% 
  ggplot(aes(x=hour, y=high_mito, colour=plotID, shape = liths))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=high_mito-se, ymax=high_mito+se), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+  
  # geom_hline(yintercept = lyso_mean, color="black", 
  #            linetype="dashed", size=1.5)+
  scale_x_continuous(name = "Hour", limits = c(-10,130), breaks=c(0,24,48,72,96,120))+  
  scale_colour_manual("plotID",
                      # labels = c("OM control","OM + coccoliths","IM control","IM + coccoliths"),
                      values = c("control im-control" = "#006699",
                                 "viral im-control" = "#FF3300",
                                 "control im-cells+liths" = "#006699",  
                                 "viral im-cells+liths" = "#FF3300"))+   
  scale_shape_manual(values = c(1,19,1,19))+
  scale_y_continuous(name = "Percentage of cells\nwith healthy mitochondria")

#source data
write.csv(IM_mito, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs4b.csv")

#need to run the photophysiology script before running the next line
ggarrange(cells, calcified, liths, fv_fm_plot, ros, sytox, align = 'hv', label.x = 0, 
          nrow = 3, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = TRUE, legend = "bottom")

ggarrange(calcified, liths, lyso, mito, align = 'hv', label.x = 0, 
          nrow = 2, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = TRUE, legend = "bottom")

#IM protliths######################################################################################################
IM_proliths <- filter(induced_media, experiment.name == "proliths_2")

IM_proliths$hour <- as.numeric(IM_proliths$hour)
IM_proliths$plotID <- as.factor(paste(IM_proliths$treatment, IM_proliths$liths, sep="-"))

IM_proliths_cells <- summarySE(IM_proliths,measurevar = "cell.concentration",  
                      groupvars = c("hour","treatment","liths","plotID"))

IM_prolith_liths <- summarySE(IM_proliths,measurevar = "lith.concentration",  
                      groupvars = c("hour","treatment","liths","plotID"))

#plotting the growth curve 
prolith_growth <- IM_proliths_cells %>% filter(treatment != "prolith supernatant") %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "untreated_liths",
                                          "proliths",
                                          "prolith supernatant"))) %>%  
  ggplot(aes(x=hour, y=cell.concentration, colour=treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=cell.concentration-se, ymax=cell.concentration+se), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144))+  
  scale_y_continuous(name = bquote(atop("Cell concentration", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("Control", "Coccoliths","Proliths","Supernatant"),
                      values = c("control" = "#FFCC00",  
                                 "untreated_liths" = "#006699",
                                 "proliths" = "#000000",
                                 "prolith supernatant" = "#CC79A7"))

#plotting the growth curve 
prolith_liths <- IM_prolith_liths %>%
  arrange(lith.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "untreated_liths",
                                                "proliths",
                                                "prolith supernatant"))) %>%  
  ggplot(aes(x=hour, y=log10(lith.concentration), colour=treatment))+  
  geom_point(size=5)+theme_bw()+ 
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2), 
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=log10(lith.concentration-se), ymax=log10(lith.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144))+  
  scale_y_continuous(name = bquote(atop("Lith concentration; log10", "liths mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("Control", "Coccoliths","Proliths","Supernatant"),
                      values = c("control" = "#FFCC00",  
                                 "untreated_liths" = "#006699",
                                 "proliths" = "#000000",
                                 "prolith supernatant" = "#CC79A7"))


#plotting the percent calcified
prolith_cal <- IM_proliths %>% filter(treatment != "prolith supernatant") %>% filter(hour < 73) %>%
  arrange(calcified) %>%
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "untreated_liths",
                                                "proliths",
                                                "prolith supernatant"))) %>%  
  ggplot(aes(x=as.factor(hour), y=calcified, colour=treatment))+  
  geom_boxplot()+
  geom_point(aes(colour =treatment), size=5, position = position_jitterdodge())+theme_bw()+ 
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  scale_x_discrete(name = "Time (h)", breaks=c(0,24,48,72,96,120,144))+  
  scale_y_continuous(name = "Percent calcified", limits =c(0,60))+
  scale_colour_manual("Treatment",
                      labels = c("Control", "Coccoliths","Proliths","Supernatant"),
                      values = c("control" = "#FFCC00",  
                                 "untreated_liths" = "#006699",
                                 "proliths" = "#000000",
                                 "prolith supernatant" = "#CC79A7"))

ggarrange(prolith_growth, prolith_liths, prolith_cal, align = 'hv', label.x = 0, 
          nrow = 2, ncol = 2, hjust = -0.35, vjust = 26, 
          common.legend = TRUE, legend = "bottom")

#The impact of boiling the induced media##########################################################################
boiled <- subset(induced_media, experiment.name == "IM_heat_killed")

boiled_cells <- boiled %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment","liths"))

boiled <- boiled_cells %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control im",  
                                          "viral im",
                                          "viral im boiled"))) %>%  
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("OM","IM","IM boiled"),
                      values = c("control im" = "#FFCC00",
                                 "viral im" = "#FF3300",
                                 "viral im boiled" = "#88CCEE"))

#IM freeze thaw cycles x3#########################################################################################
im_freeze_thaw <- subset(IM, experiment.name =="im_freeze_thaw")

im_freeze_thaw_cells <- im_freeze_thaw %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"))


freeze <- im_freeze_thaw_cells %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control im",
                                                "viral im", 
                                                "viral im freeze thaw"))) %>%  
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  theme_classic()+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("OM","IM","IM freeze thaw"),
                      values = c("control im" = "#FFCC00",
                                 "viral im" = "#FF3300",
                                 "viral im freeze thaw" = "#CC6677"))+  
  scale_shape_manual(values = c(19,19,22))

#IM protK digestion experiment###########################################################################################
im_protease <- subset(IM, experiment.name == "im_protease")

im_protease_cells <- im_protease %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"))


protease <- im_protease_cells %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control im",  
                                                "control im+protease",
                                                "viral im", 
                                                "viral im+protease"))) %>%  
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("OM", "OM +protease","IM","IM +protease"),
                      values = c("control im" = "#FFCC00",  
                                 "control im+protease" = "#FFCC00",
                                 "viral im" = "#FF3300",
                                 "viral im+protease" = "#FF3300"))+  
  scale_shape_manual(values = c(19,21,19,21))

#ion exchange columns#####################################################################################################
#this was done with 10kDa TFF concentrate that was first passed through 100kDa filter cassette
ion_screen <- subset(induced_media, experiment.name == "ion_exchange")

ion_cells <- ion_screen %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("experiment.name", "hour","treatment"), .drop = TRUE)


ion_exchange <- ion_cells %>% filter(!is.na(cell.concentration)) %>%
  filter(treatment != "control") %>%
  filter(treatment != "50kDa_concentrate") %>%
  filter(treatment != "50kDa_filtrate") %>%
  filter(treatment != "10kDa_filtrate") %>%  
  filter(treatment != "q_eluent") %>%
  filter(treatment != "s_eluent") %>%
  filter(hour < 145) %>%  
  arrange(cell.concentration) %>%  
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "10kDa_filtrate",
                                                "50kDa_concentrate",
                                                "50kDa_filtrate",
                                                "s_filtrate",  
                                                "s_eluent",  
                                                "q_filtrate",  
                                                "q_eluent"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", limits = c(-3,171), breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      values = c("control" = "#FFCC00",
                                 "s_filtrate" = "#FF3300",
                                 "s_eluent" = "#FF3300",
                                 "q_filtrate" = "#FF3300",
                                 "q_eluent" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(17,18,19,15,0))

#prepare data for export
ion_cells <- ion_cells %>% filter(!is.na(cell.concentration)) %>%
  filter(treatment != "control") %>%
  filter(treatment != "50kDa_concentrate") %>%
  filter(treatment != "50kDa_filtrate") %>%
  filter(treatment != "10kDa_filtrate") %>%  
  filter(treatment != "q_eluent") %>%
  filter(treatment != "s_eluent") %>%
  filter(hour < 145)
   
#write source data file      
# write.csv(ion_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs6c.csv")

#Plotting all of the treatments together##########################################################################
protein_character <- subset(induced_media, experiment.name == "IM_heat_killed" | experiment.name == "im_freeze_thaw" | experiment.name == "im_protease")

protein_character_cells <- protein_character %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"))

protein <- protein_character_cells %>%  
  filter(treatment != "control im+protease") %>%  
  filter(treatment != "control im") %>%  
  ggplot(aes(x= hour, y=log10(cell.concentration), colour = treatment, shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)",limits = c(-3,171), breaks=c(0,24,48,72,96,120,144,168))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+  
  scale_colour_manual(values = c("viral im" = "#FF3300",  
                                 "viral im boiled" = "#FF3300",  
                                 "viral im freeze thaw" = "#FF3300",
                                 "viral im+protease" = "#FF3300"), guide = "none")+  
  scale_shape_manual("Treatment",
                     labels = c("IM","IM boiled", "IM freeze-thaw", "IM protease"),  
                     values = c(16,15,17,18))

#source data
# write.csv(protein_character_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs6a.csv")

#IM centricon size fractionation experiment######################################################################################
im_fraction <- subset(induced_media, experiment.name == "tff_100kDa_test" | experiment.name == "im_fractionation_2" | treatment == "50kDa_concentrate" | treatment == "50kDa_filtrate")
im_fraction <- im_fraction %>% filter(cell.concentration != "NA") %>% filter(experiment.name != "hapto_screen_2")

im_fraction_cells <- im_fraction %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment","liths"))

fractionate <- im_fraction_cells %>%  
  filter(treatment != "control im") %>%
  filter(treatment != "100kDa_filtrate") %>%
  filter(treatment != "100kDa_filtrate_heated") %>%
  filter(treatment != "50kDa_concentrate") %>%
  filter(treatment != "10kDa_concentrate") %>%  
  filter(hour < 170) %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("viral im",
                                                "50kDa_concentrate",
                                                "50kDa_filtrate",
                                                "30kDa_filtrate", 
                                                "10kDa_concentrate",
                                                "10kDa_filtrate",
                                                "3kDa_filtrate"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", limits = c(-3,171), breaks=c(0,24,48,72,96,120,144,168))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("viral im" = "#FF3300",
                                 "50kDa_concentrate" = "#FF3300",
                                 "50kDa_filtrate" = "#FF3300",
                                 "30kDa_filtrate" = "#FF3300",
                                 "10kDa_concentrate" = "#FF3300",  
                                 "10kDa_filtrate" = "#FF3300",  
                                 "3kDa_filtrate" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,15,4,17,18,7,8))

#prepare data for export
im_fraction_cells <- im_fraction_cells %>%  
  filter(treatment != "control im") %>%
  filter(treatment != "100kDa_filtrate") %>%
  filter(treatment != "100kDa_filtrate_heated") %>%
  filter(treatment != "50kDa_concentrate") %>%
  filter(treatment != "10kDa_concentrate") %>%
  filter(hour < 170)
#write source data file
write.csv(im_fraction_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs6b.csv")

ggarrange(protein, fractionate, ion_exchange, align = 'hv', label.x = 0, 
          nrow = 2, ncol = 2)

#Does IM kill other haptophytes###################################################################################
hapto_screen <- subset(induced_media, experiment.name == "hapto_screen")
hapto_screen <- subset(hapto_screen, strain_id == "g_ocean")

hapto_screen_2 <- subset(induced_media, experiment.name == "hapto_screen_2")

hapto_screen_2 <- rbind(hapto_screen, hapto_screen_2)

hapto_screen_2$plotID <- as.factor(paste(hapto_screen_2$strain_id, hapto_screen_2$treatment, sep="-"))

hapto_screen_cells <- hapto_screen_2 %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("hour","strain_id","treatment","plotID"))


hapto_screen_cells %>% 
  filter(treatment != "50kDa_filtrate") %>%
  filter(hour != "24") %>%
  filter(hour < 121) %>%
  arrange(cell.concentration) %>%
  mutate(plotID = factor(plotID, levels=c("374-control im",
                                                "374-viral im",
                                                "iso-viral im",
                                                "tiso-viral im",
                                                "phaeo-viral im",
                                          "g_ocean-viral im"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=plotID, shape = plotID))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+
  scale_x_continuous(name = "Hour", limits = c(0,130), breaks=c(0,24,48,72,96,120))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("CCMP374 Control","CCMP374 Viral IM","Isochrysis Viral IM", "Tisochrysis Viral IM", "Phaeocystis Viral IM", "Gephyrocapsa Viral IM"),
                      values = c("374-control im" = "#006699",
                                 "374-viral im" = "#FF3300",
                                 "iso-viral im" = "#FF3300",
                                 "tiso-viral im" = "#FF3300",
                                 "phaeo-viral im" = "#FF3300",
                                 "g_ocean-viral im" = "#FF3300"))+  
  scale_shape_manual(values = c(19,19,15,17,18,3))

#prep for export
hapto_screen_cells <- hapto_screen_cells %>% 
  filter(treatment != "50kDa_filtrate") %>%
  filter(hour != "24") %>%
  filter(hour < 121)

#source data
# write.csv(hapto_screen_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs5.csv")

#Testing the dynamics after putting the large culture through the TFF pipeline and ammonium sulfate fractionation############################### 
#using a 100kDa membrane to remove the virus and a 10kDa to concentrate the protein
tff_test <- subset(induced_media, experiment.name == "tff_100kDa_test")

tff_cells <- tff_test %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("hour","treatment"), .drop = TRUE)

tff_protein_concentrations <- filter(tff_protein_concentrations, sample != "std")
boxplot(concentration~sample, data = tff_protein_concentrations, ylab = "Protein concentration")

ammonium_screen <- subset(induced_media, experiment.name == "ammonium_screen_2")

ammonium_screen$precipitation_media <- rep(c("fsw","tris"), each = 240, length=480)

ammonium_cells <- ammonium_screen %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("experiment.name", "hour","treatment"), .drop = TRUE)


fsw <- ammonium_cells %>% filter(!is.na(cell.concentration)) %>%
  arrange(cell.concentration) %>%  
  filter(treatment != "control") %>%
  # filter(treatment != "20_ammonium") %>%
  # filter(treatment != "100_ammonium") %>%
  filter(experiment.name == "ammonium_screen_2") %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "10kDa_filtrate",
                                                "10kDa_concentrate",
                                                "20_ammonium",  
                                                "40_ammonium",  
                                                "60_ammonium",  
                                                "80_ammonium",  
                                                "100_ammonium"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168,192,216))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",  
                      values = c("10kDa_filtrate" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",  
                                 "20_ammonium" = "#FF3300",  
                                 "40_ammonium" = "#FF3300",  
                                 "60_ammonium" = "#FF3300",  
                                 "80_ammonium" = "#FF3300",  
                                 "100_ammonium" = "#FF3300"), guide = "none")+  
  scale_shape_manual("Treatment", values = c(19,0,1,2,3,4,5))

#prepare data for export
ammonium_cells <- ammonium_cells %>% filter(!is.na(cell.concentration)) %>%
  filter(treatment != "control") %>%
  filter(experiment.name == "ammonium_screen_2")

#write source data file
# write.csv(ammonium_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs10.csv")

#Do coccoliths increase encounters of vesicles and cells##############################################################################################################################
#vesicle concentrations in response to induced media
#for vesicle concentrations combined from both 100 and 50kDa fitlrates
vesicle_measurements <- filter(induced_media, !is.na(vesicle.concentration))
#only using the control and 10kDa treatments
vesicle_measurements <- vesicle_measurements %>% subset(treatment == "control" | treatment == "10kDa_concentrate")
#summarizing the vesicle concentrations
vesicle_measurements <- vesicle_measurements %>% filter(experiment.name != "proliths_vesicle")
vesicle_measurements <- vesicle_measurements %>% filter(experiment.name != "ammonium_screen")

vesicle_measurements_sum <- vesicle_measurements %>%
  summarySE(measurevar = "vesicle.concentration.expanded",  
            groupvars = c("hour","treatment"), na.rm= TRUE, .drop = TRUE)
#plotting
vesicle_plot <- vesicle_measurements_sum %>%
  filter(hour < 192) %>%
  ggplot(aes(x=hour, y=log10(vesicle.concentration.expanded), colour=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(vesicle.concentration.expanded-se), ymax=log10(vesicle.concentration.expanded+se)), width= 5, size=1)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)",  
                     breaks=c(0,24,48,72,96,120,144,168,192,216,240,264,288,312,336,360,384))+  
  scale_y_continuous(name = bquote(atop("Vesicle concentration; log10", "vesicles mL"^-1)))+
  scale_colour_manual("Treatment", 
                      values = c("control" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300"), guide = "none")

#prep file
vesicle_measurements_sum <- vesicle_measurements_sum %>%
  filter(hour < 192)

#source data
# write.csv(vesicle_measurements_sum, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig5b.csv")

#note that the spike at 72 h is driven by a burst of smaller particles in one of the experiments, these are smaller than the typical vesicles that I've been measuring, at 96 you can see there is a transition from the smaller particles
#for the cell concentrations as requested by the reviewer
cell_measurements_sum <- vesicle_measurements %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("hour","treatment"), .drop = TRUE)

experiment_fig5_panel_c_cell_con_sum <- experiment_fig5_panel_c %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("hour","treatment"), .drop = TRUE)
#ploting
cell_plot <- cell_measurements_sum %>%
  filter(hour < 192) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=1)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)",  
                     breaks=c(0,24,48,72,96,120,144,168,192,216,240,264,288,312,336,360,384))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment", 
                      values = c("control" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300"), guide = "none")

#prep file
cell_measurements_sum <- cell_measurements_sum %>%
  filter(hour < 192)

#source data
# write.csv(cell_measurements_sum, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig5a.csv")

ggarrange(cell_plot, vesicle_plot, nrow = 2, ncol = 2, align = "h")

#Vesicles were harvested from cells and incubated for 2 h with free coccoliths, they were then added to cells at a ratio of 1 vesicle for every 10 cells
#a control with no coccoliths was also set up to compare dynamics
#I was also interested in what component was more important in initiating cell death (protein or vesicles) so I also have two treatments (with and without coccoliths)
#that have been spiked with the same volume I added for the vesicle treatment except I removeed the vesicles (using 0.02um anotop) leaving just the dissolved protein
vesicle_liths <- subset(induced_media, experiment.name == "vesicle_liths")
vesicle_liths$plotID <- as.factor(paste(vesicle_liths$treatment, vesicle_liths$liths, sep="-"))

vl_cells <- vesicle_liths %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("experiment.name", "hour","treatment",  "liths","plotID"), .drop = TRUE)


vesicle_removal <- vl_cells %>% filter(!is.na(cell.concentration)) %>%
  arrange(cell.concentration) %>%  
  filter(liths != "cells+liths") %>%
  mutate(plotID = factor(plotID, levels=c("control-control",  
                                                "control-cells+liths",
                                                "ev_filtrate-control",
                                                "ev_filtrate-cells+liths",
                                                "evs-control",  
                                                "evs-cells+liths"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=plotID, shape = plotID))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+
  scale_x_continuous(name = "Time (h)",  
                     breaks=c(0,24,48,72,96,120,144,168,192,216,240,264,288,312,336,360,384))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
scale_colour_manual("Treatment", 
                    values = c("control-control" = "#FFCC00",
                               "ev_filtrate-control" = "#D55E00",
                               "evs-control" = "#D55E00"), guide = "none")+
scale_shape_manual("Treatment", labels = c("Control","Vesicle-free Media","Vesicles"), values = c(19,17,2))

#prep file
vl_cells <- vl_cells %>% filter(liths != "cells+liths")

#source data
# write.csv(vl_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig5f.csv")

ggarrange(vesicle_plot, vesicle_removal, nrow = 2, ncol = 2, align = "h")

#Protease inhibitor experiments############################################################################################################################
#these experiments were based on the finding of a serine protease in all of the ammonium sulfate fractions, I hypothesized that using a serine protease
#inhibitor would knock out the activity of the protein in the IM, I used E-64 which is a cysteine protease inhibitor and ABESF which is a serine protease inhibitor
protease_inhibitors <- induced_media %>% subset(experiment.name == "protease_inhibitors_pre_inc" |  
                                                  experiment.name == "protease_inhibitors_vesicles" |  
                                                  experiment.name == "protease_inhibitors_caspase" |  
                                                  experiment.name == "protease_inhibitors_everything")

protease_inhibitors_cells <- protease_inhibitors %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"))

protease_cells <- protease_inhibitors_cells %>%  
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "10kDa_concentrate",
                                                "e_64_10kDa_concentrate",
                                                "aebsf_10kDa_concentrate"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)), breaks = c(3,4,4,5,6), limits = c(3,6.5))+
  scale_colour_manual("Treatment",
                      # labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("control" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "e_64_10kDa_concentrate" = "#FF3300",
                                 "aebsf_10kDa_concentrate" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,17,3,5))

#source data
# write.csv(protease_inhibitors_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig6a.csv")

#filter 120 h timepoint
timepoint_120 <- protease_inhibitors %>% filter(hour == "120")
timepoint_120 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$cell.concentration)))#normal
bartlett.test(log10(cell.concentration)~interaction(treatment,experiment.name), data=timepoint_120)#variances of group are not really different after log transformation
#boxplots
boxplot(cell.concentration ~ treatment, data = timepoint_120, xlab = "treatment", ylab = "cell concentration")
#anova
mod1 <- aov(log10(cell.concentration) ~ experiment.name * treatment, data = timepoint_120)
summary(mod1)
hist(residuals(mod1), col = "darkgrey")
plot(fitted(mod1), residuals(mod1))
TukeyHSD(mod1, which = "treatment")
HSD.test(mod1, trt = c("treatment"), console = TRUE)

#filter 144 h timepoint
timepoint_144 <- protease_inhibitors %>% filter(hour == "144")
timepoint_144 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$cell.concentration)))#normal
bartlett.test(log10(cell.concentration)~interaction(treatment,experiment.name), data=timepoint_144)#variances of group are not really different after log transformation

#boxplots
boxplot(cell.concentration ~ treatment, data = timepoint_144, xlab = "treatment", ylab = "cell concentration")
#anova
mod2 <- aov(log10(cell.concentration) ~ experiment.name + treatment, data = timepoint_144)
summary(mod2)
hist(residuals(mod2), col = "darkgrey")
plot(fitted(mod2), residuals(mod2))
TukeyHSD(mod2, which = "treatment")
HSD.test(mod2, trt = c("treatment"), console = TRUE)

#filtering 169 h timepoint
timepoint_168 <- protease_inhibitors %>% filter(hour == "168")
timepoint_168 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$cell.concentration)))#normal
bartlett.test(log10(cell.concentration)~interaction(treatment,experiment.name), data=timepoint_168)#variances of group are not really different after log transformation
#boxplots
boxplot(cell.concentration ~ treatment, data = timepoint_168, xlab = "treatment", ylab = "cell concentration")
#anova
mod3 <- aov(log10(cell.concentration) ~ experiment.name + treatment, data = timepoint_168)
summary(mod3)
hist(residuals(mod3), col = "darkgrey")
plot(fitted(mod3), residuals(mod3))
TukeyHSD(mod3, which = "treatment")
HSD.test(mod3, trt = c("treatment"), console = TRUE)

#capsase data
protease_inhibitors_caspase <- induced_media %>% subset(experiment.name == "protease_inhibitors_caspase" | experiment.name == "protease_inhibitors_everything")

#removes the rows with no data otherwise the summarysE function does not work properly
protease_inhibitors_caspase <- protease_inhibitors_caspase[!is.na(protease_inhibitors_caspase$caspase_perc_pos), ]

protease_inhibitors_capase_sum <- protease_inhibitors_caspase %>% summarySE(measurevar = "caspase_perc_pos", groupvars = c("hour","treatment"))

protease_caspase <- protease_inhibitors_capase_sum %>%  
  arrange(caspase_perc_pos) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "10kDa_concentrate",
                                                "e_64_10kDa_concentrate",
                                                "aebsf_10kDa_concentrate"))) %>%
  ggplot(aes(x=hour, y=caspase_perc_pos, colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=caspase_perc_pos-se, ymax=caspase_perc_pos+se), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = "Caspase positive cells")+
  scale_colour_manual("Treatment",
                      # labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("control" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "e_64_10kDa_concentrate" = "#FF3300",
                                 "aebsf_10kDa_concentrate" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,17,3,5))

#source data
# write.csv(protease_inhibitors_capase_sum, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig6c.csv")

#stats for the caspase measurements at each of the critical timepoints 120-168 h
timepoint_144 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$caspase_perc_pos)))#normal
bartlett.test(log10(caspase_perc_pos)~treatment, data=timepoint_120)#variances of group are not really different after log transformation

mod1 <- aov(log10(caspase_perc_pos) ~ treatment, data = timepoint_120)
summary(mod1)
hist(residuals(mod1), col = "darkgrey")
plot(fitted(mod1), residuals(mod1))
TukeyHSD(mod1, which = "treatment")
HSD.test(mod1, trt = c("treatment"), console = TRUE)

bartlett.test(caspase_perc_pos~interaction(treatment,experiment.name), data=timepoint_144)#variances of group are not really different after log transformation
mod2 <- aov(caspase_perc_pos ~ experiment.name + treatment, data = timepoint_144)
summary(mod2)
hist(residuals(mod2), col = "darkgrey")
plot(fitted(mod2), residuals(mod2))
TukeyHSD(mod2, which = "treatment")
HSD.test(mod2, trt = c("treatment"), console = TRUE)

bartlett.test(caspase_perc_pos~interaction(treatment,experiment.name), data=timepoint_168)#variances of group are not really different after log transformation
mod3 <- aov(caspase_perc_pos ~ experiment.name * treatment, data = timepoint_168)
summary(mod3)
hist(residuals(mod3), col = "darkgrey")
plot(fitted(mod3), residuals(mod3))
TukeyHSD(mod3, which = "treatment")
HSD.test(mod3, trt = c("treatment"), console = TRUE)

#looks at the impact on vesicle production
protease_inhibitors_vesicles <- induced_media %>% subset(experiment.name == "protease_inhibitors_vesicles" | experiment.name == "protease_inhibitors_caspase" | experiment.name == "protease_inhibitors_everything")

protease_inhibitors_vesicles <- protease_inhibitors_vesicles[!is.na(protease_inhibitors_vesicles$vesicle.concentration), ]


protease_inhibitors_vesicles_sum <- protease_inhibitors_vesicles %>% summarySE(measurevar = "vesicle.concentration.expanded", groupvars = c("hour","treatment"), na.rm = TRUE)

protease_vesicles <- protease_inhibitors_vesicles_sum %>%  
  arrange(vesicle.concentration.expanded) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "10kDa_concentrate",
                                                "e_64_10kDa_concentrate",
                                                "aebsf_10kDa_concentrate"))) %>%
  ggplot(aes(x=hour, y=log10(vesicle.concentration.expanded), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(vesicle.concentration.expanded-se), ymax=log10(vesicle.concentration.expanded+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Vesicle concentration; log10", "vesicles mL"^-1)), breaks = c(6,6.5,7,7.5,8,8.5), limits = c(6,8.5))+
  scale_colour_manual("Treatment",
                      # labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("control" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "e_64_10kDa_concentrate" = "#FF3300",
                                 "aebsf_10kDa_concentrate" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,17,3,5))

#source data
# write.csv(protease_inhibitors_vesicles_sum, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig6b.csv")

#stats for vesicle measurements for the same time points, 120-168h
timepoint_120 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$vesicle.concentration)))#normal
bartlett.test(log10(vesicle.concentration)~treatment, data=timepoint_120)#variances of group are not really different after log transformation

mod1 <- aov(log10(vesicle.concentration.expanded) ~ treatment, data = timepoint_120)
summary(mod1)
plot(mod1)
hist(residuals(mod1), col = "darkgrey")
plot(fitted(mod1), residuals(mod1))
TukeyHSD(mod1, which = "treatment")
HSD.test(mod1, trt = c("treatment"), console = TRUE)

#for the 144 h time vesicle time point
timepoint_144 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$vesicle.concentration)))#normal
bartlett.test(log10(vesicle.concentration)~treatment, data=timepoint_144)#variances of group are not really different after log transformation

welch_anova_test(log10(vesicle.concentration.expanded)~as.factor(treatment), data=timepoint_144)
timepoint_144$log10_vesicle_con_expanded <- log10(timepoint_144$vesicle.concentration.expanded)
games_howell_test(log10_vesicle_con_expanded~treatment, data=timepoint_144)

#filtering 168 h timepoint
timepoint_168_vesicle <- timepoint_168 %>% filter(experiment.name == "protease_inhibitors_caspase" | experiment.name == "protease_inhibitors_everything")
timepoint_168_vesicle$log10_vesicle_concentrations <- as.numeric(log10(timepoint_168_vesicle$vesicle.concentration))
timepoint_168_vesicle$log10_vesicle_concentrations <- as.numeric(timepoint_168_vesicle$log10_vesicle_concentrations)

timepoint_168 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$log10_vesicle_concentrations)))#normal
bartlett.test(log10(vesicle.concentration)~interaction(treatment,experiment.name), data=timepoint_168)#variances of group are not really different after log transformation
bartlett.test(log10(vesicle.concentration.expanded)~treatment, data=timepoint_168)#variances of group are not really different after log transformation
boxplot(log10(vesicle.concentration.expanded)~treatment, data=timepoint_168)

kruskal.test(log10(vesicle.concentration.expanded) ~ treatment, data = timepoint_168)
dunn_res <- dunnTest(log10(vesicle.concentration.expanded) ~ treatment,
         data=timepoint_168,
         method="bh")    

cldList(P.adj ~ Comparison,
  data = dunn_res$res,
  threshold = 0.05)

#need to run the last section before this
ggarrange(protease_cells, protease_caspase, protease_vesicles, align = 'hv', nrow = 2, ncol = 2, common.legend = FALSE)

#I'm going to add the "subtilisin" inhibitor in this section too- aka PF-429242
#this was the first experiment I did with this inhibitor using 20uM and 10uM spikes at t0 and 72h
sub_inhibitor_exp <- induced_media %>% filter(experiment.name == "pf429242_inhibitor")

sub_inhibitor_exp_cells <- sub_inhibitor_exp %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"))

sub_inhibitor_exp_cells %>%  
  # filter(hour != "120") %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "pf429242",
                                                "10kDa_concentrate",
                                                "10kDa_concentrate_pf429242"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      values = c("control" = "#FFCC00",  
                                 "pf429242" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "10kDa_concentrate_pf429242" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,1,17,2,6,7,8))

timepoint_168_experiment1 <- sub_inhibitor_exp %>% filter(hour == "168")
timepoint_168_experiment1 %>% group_by(treatment) %>% do(tidy(shapiro.test(.$cell.concentration)))#normal
bartlett.test(log10(cell.concentration)~treatment, data=timepoint_168_experiment1)#variances of group are not really different after log transformation

#anova
mod1 <- aov(log10(cell.concentration) ~ treatment, data = timepoint_168_experiment1)
summary(mod1)
hist(residuals(mod1), col = "darkgrey")
plot(fitted(mod1), residuals(mod1))
TukeyHSD(mod1, which = "treatment")
HSD.test(mod1, trt = c("treatment"), console = TRUE)

#the second and third experiments performed with the PF429242
sub_inhibitor_exp_2_3_4_5 <- induced_media %>% subset(experiment.name == "pf429242_inhibitor_2" | experiment.name == "pf429242_inhibitor_3" | experiment.name == "pf429242_inhibitor_4" | experiment.name == "pf429242_inhibitor_5")

sub_inhibitor_exp_2_3_4_5_cells <- sub_inhibitor_exp_2_3_4_5 %>%  summarySE(measurevar = "cell.concentration", groupvars = c("experiment.name", "hour","treatment"))

sub_inhibitor_exp_2_3_4_5_bodipy <- sub_inhibitor_exp_2_3_4_5 %>%  summarySE(measurevar = "percent_bodipy", groupvars = c("experiment.name", "hour","treatment"))

sub_inhibitor_exp_2_3_4_5_vesicles <- sub_inhibitor_exp_2_3_4_5 %>%  summarySE(measurevar = "vesicle.concentration.expanded", groupvars = c("hour","treatment"), na.rm = TRUE)


sub_inhibitor_exp_2_3_4_5_cells %>%  
  # filter(hour != "120") %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "pf429242",
                                                "10kDa_concentrate",
                                                "10kDa_concentrate_pf429242"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      values = c("control" = "#FFCC00",  
                                 "pf429242" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "10kDa_concentrate_pf429242" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,1,17,2,6,7,8))+
  facet_wrap(~experiment.name, dir ="v")

#soure data
write.csv(sub_inhibitor_exp_2_3_4_5_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs12c-d.csv")

#plots the last experiment with the subtilisin inhibitor
pf_cells <- sub_inhibitor_exp_2_3_4_5_cells %>%  
  filter(experiment.name == "pf429242_inhibitor_5") %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "pf429242",
                                                "10kDa_concentrate",
                                                "10kDa_concentrate_pf429242"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)), breaks = c(3,4,4,5,6), limits = c(3,6.5))+
  scale_colour_manual("Treatment",
                      values = c("control" = "#FFCC00",  
                                 "pf429242" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "10kDa_concentrate_pf429242" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,1,17,2,6,7,8))

#prep file
sub_inhibitor_exp5 <- sub_inhibitor_exp_2_3_4_5_cells %>%  
  filter(experiment.name == "pf429242_inhibitor_5")

#source data
# write.csv(sub_inhibitor_exp5, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig6d.csv")

#stats on the cells subsetting the different timepoints
timepoint_120_subtilisin <- sub_inhibitor_exp_2_3_4_5 %>% filter(experiment.name == "pf429242_inhibitor_5" & hour == "120")
timepoint_144_subtilisin <- sub_inhibitor_exp_2_3_4_5 %>% filter(experiment.name == "pf429242_inhibitor_5" & hour == "144")
timepoint_168_subtilisin <- sub_inhibitor_exp_2_3_4_5 %>% filter(experiment.name == "pf429242_inhibitor_5" & hour == "168")

#anova for cell concentration looks good
timepoint_168_subtilisin %>% group_by(treatment) %>% do(tidy(shapiro.test(.$cell.concentration)))#normal
bartlett.test(log10(cell.concentration)~treatment, data=timepoint_168_subtilisin)#variances of group are not really different after log transformation

mod1 <- aov(log10(cell.concentration) ~ treatment, data = timepoint_168_subtilisin)
summary(mod1)
hist(residuals(mod1), col = "darkgrey")
plot(fitted(mod1), residuals(mod1))
TukeyHSD(mod1, which = "treatment")
HSD.test(mod1, trt = c("treatment"), console = TRUE)

timepoint_168_subtilisin %>% group_by(treatment) %>% do(tidy(shapiro.test(.$vesicle.concentration)))#normal
bartlett.test(log10(vesicle.concentration)~treatment, data=timepoint_168_subtilisin)#variances of group are not really different after log transformation

mod2 <- aov(log10(vesicle.concentration) ~ treatment, data = timepoint_168_subtilisin)
summary(mod2)
hist(residuals(mod2), col = "darkgrey")
plot(fitted(mod2), residuals(mod2))
TukeyHSD(mod2, which = "treatment")
HSD.test(mod2, trt = c("treatment"), console = TRUE)

pf_bodipy <- sub_inhibitor_exp_2_3_4_5_bodipy %>%  
  filter(experiment.name == "pf429242_inhibitor_5") %>%
  arrange(percent_bodipy) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "pf429242",
                                                "10kDa_concentrate",
                                                "10kDa_concentrate_pf429242"))) %>%
  ggplot(aes(x=hour, y=percent_bodipy, colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=percent_bodipy-se, ymax=percent_bodipy+se), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = "Percent Bodipy")+
  scale_colour_manual("Treatment",
                      values = c("control" = "#FFCC00",  
                                 "pf429242" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "10kDa_concentrate_pf429242" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,1,17,2,6,7,8))

#prep file
sub_inhibitor_exp5_bodipy <- sub_inhibitor_exp_2_3_4_5_bodipy %>%  
  filter(experiment.name == "pf429242_inhibitor_5")


#source data
write.csv(sub_inhibitor_exp5_bodipy, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig6f.csv")

#stats
timepoint_168_subtilisin %>% group_by(treatment) %>% do(tidy(shapiro.test(.$percent_bodipy)))#normal
bartlett.test(percent_bodipy~treatment, data=timepoint_168_subtilisin)#variances of group are not really different after log transformation

welch_anova_test(percent_bodipy~as.factor(treatment), data=timepoint_168_subtilisin)
ght <- games_howell_test(percent_bodipy~treatment, data=timepoint_168_subtilisin)
p_values <- ght$p.adj
names(p_values) <- paste0(ght$group1, "-", ght$group2)

cld <- multcompLetters(p_values, threshold = 0.05)


pf_vesicles <- sub_inhibitor_exp_2_3_4_5_vesicles %>%  
  arrange(vesicle.concentration.expanded) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "pf429242",
                                                "10kDa_concentrate",
                                                "10kDa_concentrate_pf429242"))) %>%
  ggplot(aes(x=hour, y=log10(vesicle.concentration.expanded), colour=treatment, shape=treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),
        axis.text.x = element_text(colour = "black", size = 16),
        axis.text.y = element_text(colour = "black", size = 16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.2, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  geom_errorbar(aes(ymin=log10(vesicle.concentration.expanded-se), ymax=log10(vesicle.concentration.expanded+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Vesicle concentration; log10", "cells mL"^-1)), breaks = c(6,6.5,7,7.5,8,8.5), limits = c(6,8.5))+
  scale_colour_manual("Treatment",
                      values = c("control" = "#FFCC00",  
                                 "pf429242" = "#FFCC00",
                                 "10kDa_concentrate" = "#FF3300",
                                 "10kDa_concentrate_pf429242" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,1,17,2,6,7,8))

#source data
# write.csv(sub_inhibitor_exp_2_3_4_5_vesicles, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig6e.csv")

#anova for vesicle concentration looks good
bartlett.test(log10(vesicle.concentration.expanded)~treatment, data=timepoint_168_subtilisin)#variances of group are not really different after log transformation
sub_ves_anova <- aov(log10(vesicle.concentration.expanded)~treatment, data = timepoint_168_subtilisin)
summary(sub_ves_anova)
hist(residuals(sub_ves_anova), col = "darkgrey")
plot(fitted(mod1), residuals(mod1))
TukeyHSD(sub_ves_anova, which = "treatment")
HSD.test(sub_ves_anova, trt = c("treatment"), console = TRUE)

ggarrange(protease_cells, pf_cells, protease_vesicles, pf_vesicles, protease_caspase, pf_bodipy, nrow = 3, ncol = 2, align = "hv", common.legend = FALSE)

#here is the toxicity assay to test the concentrations to use
cell_viability <- cell_viability_anti_viral %>% filter(date == "10/14/22")
cell_viability <- data.table(cell_viability)
setkey(cell_viability, "row")
cell_viability[, ratio_loss:= (od_blank_sub/(od_blank_sub[match("0", pf429242_concentration)])), by= c("row")]

cell_viability$percent_loss <- cell_viability$ratio_loss * 100

cell_viability_percent_sum <- cell_viability %>%  summarySE(measurevar = "percent_loss", groupvars = c("pf429242_concentration"))

cell_viability_percent_sum$log_concentration <- log10(cell_viability_percent_sum$pf429242_concentration)


cell_viability_percent_sum %>% filter(pf429242_concentration < 126) %>%
  ggplot(aes(x=pf429242_concentration, y=percent_loss))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=percent_loss-se, ymax=percent_loss+se), width= 5, size=0.75)+  
  stat_smooth(method = 'nls', formula = 'y~a*exp(b*x)',
              method.args = list(start=c(a=110, b=0)), se=FALSE)

anti_viral <- cell_viability_anti_viral %>% filter(date =="10/24/22")
anti_viral <- data.table(anti_viral)
setkey(anti_viral, "row")
anti_viral[, ratio_loss:= (od_blank_sub/(od_blank_sub[match("0", pf429242_concentration)])), by= c("row")]

anti_viral$percent_loss <- anti_viral$ratio_loss * 100

anti_viral_percent_sum <- anti_viral %>%  summarySE(measurevar = "percent_loss", groupvars = c("pf429242_concentration"))

anti_viral_percent_sum$log_concentration <- log10(anti_viral_percent_sum$pf429242_concentration)

anti_viral_percent_sum %>% 
  # filter(pf429242_concentration < 126) %>%
  ggplot(aes(x=pf429242_concentration, y=percent_loss))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=percent_loss-se, ymax=percent_loss+se), width= 5, size=0.75)+  
  stat_smooth(method = 'nls', formula = 'y~a*exp(b*x)',
              method.args = list(start=c(a=110, b=0)), se=FALSE)

#using the "drda" package to calculate the cc50 of pf429242, they have a few models that I can run, and then I will check the fit
#logistic5, logistic4, logistic2,
#gompertz, loglogistic5, loglogistic4, loglogistic2, or loggompertz

model_5 <- drda(cell_viability$percent_loss ~ cell_viability$pf429242_concentration, mean_function = "loglogistic5")
model_6 <- drda(cell_viability$percent_loss ~ cell_viability$pf429242_concentration, mean_function = "loglogistic4")
model_8 <- drda(cell_viability$percent_loss ~ cell_viability$pf429242_concentration, mean_function = "loggompertz")

anova(model_6)

summary(model_6)

effective_dose(model_6, 90, type = "absolute", level = 0.95)

cell_viability_percent_sum %>% filter(pf429242_concentration > 0) 
drda(percent_loss ~ log_concentration)

#plot the dose response curve
dose_response_curve <- as_ggplot(plot(model_6, base = "10",xlim = c(0.1,1000), xlab = "PF-429242 concentration (uM)", ylab = "Percent loss"))

#Do I see reduced activity with the prolith filtrate?#####################################################
proliths_exp2 <- induced_media %>% subset(experiment.name == "proliths_2") 

proliths_exp2_cells <- proliths_exp2 %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"))
proliths_exp2_liths <- proliths_exp2 %>%  summarySE(measurevar = "calcified", groupvars = c("hour","treatment"))

prolith_1 <- proliths_exp2_cells %>%  
  filter(treatment != "prolith supernatant") %>%
  # filter(plotID != "proliths-cells+liths") %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "untreated_liths",
                                                "proliths"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape =treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", limits = c(0,168), breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      # labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("control" = "#FFCC00",
                                 "untreated_liths" = "#006699",
                                 "proliths" = "#56B4E9"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,19,19))

#prep data
proliths_exp2_cells <- proliths_exp2_cells %>%  
  filter(treatment != "prolith supernatant")

#source data
# write.csv(proliths_exp2_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs9a.csv")

prolith_2 <- proliths_exp2_liths %>%
  filter(hour < 120) %>%
  filter(treatment != "prolith supernatant") %>%
  # filter(plotID != "proliths-cells+liths") %>%
  arrange(calcified) %>%
  mutate(treatment = factor(treatment, levels=c("control",
                                                "untreated_liths",
                                                "proliths"))) %>%
  ggplot(aes(x=hour, y=calcified, colour=treatment, shape =treatment))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=calcified-se, ymax=calcified+se), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", limits = c(0,168), breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = "Percent calcified")+
  scale_colour_manual("Treatment",
                      # labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("control" = "#FFCC00",
                                 "untreated_liths" = "#006699",
                                 "proliths" = "#56B4E9"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,19,19))

#prep data
proliths_exp2_liths <- proliths_exp2_liths %>%
  filter(hour < 120) %>%
  filter(treatment != "prolith supernatant")

#source data
# write.csv(proliths_exp2_liths, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs9b.csv")

proliths_exp3 <- induced_media %>% subset(experiment.name == "proliths_exp3")

proliths_exp3$plotID <- as.factor(paste(proliths_exp3$treatment, proliths_exp3$liths, sep="-"))


proliths_exp3_cells <- proliths_exp3 %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment", "plotID"))

prolith_3 <- proliths_exp3_cells %>%  
  filter(hour != "24") %>%
  filter(hour < 169) %>%
  filter(plotID != "control-cells+liths") %>%
  filter(plotID != "prolith_supernatant-supernatant") %>%
  filter(plotID != "proliths-cells+liths") %>%
  arrange(cell.concentration) %>%
  mutate(plotID = factor(plotID, levels=c("control-control",
                                                "10kDa_concentrate-control",
                                                "prolith_filtrate-control"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=plotID, shape=plotID))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      # labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("control-control" = "#FFCC00",
                                 "10kDa_concentrate-control" = "#FF3300",
                                 "prolith_filtrate-control" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(19,19,15))

#prep data
proliths_exp3_cells <- proliths_exp3_cells %>%  
  filter(hour != "24") %>%
  filter(hour < 169) %>%
  filter(plotID != "control-cells+liths") %>%
  filter(plotID != "prolith_supernatant-supernatant") %>%
  filter(plotID != "proliths-cells+liths")

#source data
# write.csv(proliths_exp3_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs9c.csv")

proliths_exp2_liths

prolith_exps <- induced_media %>% subset(experiment.name == "proliths_exp4" | experiment.name == "proliths_exp5")

prolith_exps$plotID <- as.factor(paste(prolith_exps$percent_im_added, prolith_exps$treatment, sep="-"))

prolith_exps_cells <- prolith_exps %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment", "liths", "plotID"), na.rm = TRUE)

prolith_exps_cells$plotID_2 <- as.factor(paste(prolith_exps_cells$plotID, prolith_exps_cells$liths, sep="-"))


prolith_4 <- prolith_exps_cells %>%  
  filter(hour != "144") %>%  
  filter(hour != "216") %>%  
  filter(hour != "192") %>%
  filter(plotID != "0-control") %>%
  filter(plotID != "0.5-10kDa_concentrate") %>%
  arrange(cell.concentration) %>%
  # mutate(treatment = factor(treatment, levels=c("control",
  #                                               "10kDa_concentrate",
  #                                               "e_64_10kDa_concentrate",
  #                                               "aebsf_10kDa_concentrate"))) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=plotID, shape = plotID_2))+  
  geom_point(size=5)+theme_bw()+    
  theme_classic()+  
  theme(legend.position = c(0.2,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=plotID,group=plotID_2))+ 
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192,216,240))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      # labels = c("Viral IM", "50kDa Concentrate","50kDa Filtrate", "30kDa Filtrate", "10kDa Concentrate", "10kDa Filtrate", "3kDa Filtrate"),
                      values = c("2-10kDa_concentrate" = "#FF3300",
                                 "8-10kDa_concentrate" = "#FF3300"), guide = "none")+
  scale_shape_manual("Treatment", values = c(1,19,2,17))

#prep data
prolith_exps_cells <- prolith_exps_cells %>%  
  filter(hour != "144") %>%  
  filter(hour != "216") %>%  
  filter(hour != "192") %>%
  filter(plotID != "0-control") %>%
  filter(plotID != "0.5-10kDa_concentrate")

#source data
# write.csv(prolith_exps_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs9d.csv")


ggarrange(prolith_1, prolith_2, prolith_3, prolith_4, nrow = 2, ncol = 2, align = 'hv')

#subtilisin expression in NAAMES field data#####################################################################################
#these are averages from the original host:virus tpm ratios
subtilisin_exp_updated <- subtilisin_exp_updated %>% mutate(ratiovirushost.1 = ifelse(day >= 1 & day <= 2, "0.002284014",  
                                                                      ifelse(day >= 3 & day <= 3.9, "0.001549869",  
                                                                             ifelse(day >= 4 & day <= 6, "0.00065582",
                                                                                    ifelse(day >= 7 & day <= 9, "0.000116554","")))))
subtilisin_exp_updated$ratiovirushost <- as.numeric(subtilisin_exp_updated$ratiovirushost)
subtilisin_exp_updated <- subtilisin_exp_updated %>% mutate(cast_type = ifelse(day >= 1 & day <= 2, "biology",  
                                                                                      ifelse(day >= 3 & day <= 3.9, "biology",  
                                                                                             ifelse(day >= 4 & day <= 6, "incubation",
                                                                                                    ifelse(day >= 7 & day <= 9, "incubation","")))))

#what does the ratio of virus:host look like as we occupied station 4
boxplot(ratiovirushost~day, data = subtilisin_exp)
virus_host_stats <- subset(subtilisin_exp, day == "1" | day == "7")
boxplot(ratiovirushost~day, data = virus_host_stats)
#testing for normality
day1 <- subset(virus_host_stats, day == "1")
shapiro.test(day1$ratiovirushost)#normal
day7 <- subset(virus_host_stats, day =="7")
shapiro.test(day7$ratiovirushost)#normal
#t.test for significance
t.test(virus_host_stats$ratiovirushost,virus_host_stats$day, paired = TRUE)

#filtering the significantly expressed serine proteases
subtilisin_exp_updated %>% filter(pvalue < 0.049) %>% ggplot(aes(x= ratiovirushost.1, y = log2FoldChange))+geom_point(size = 5)

subtilisin_exp_sig <- subtilisin_exp_updated %>% filter(pvalue < 0.05)

virus_host <- subtilisin_exp %>% subset(day == "1" | day == "7") %>% ggplot(aes(x=as.factor(day), y=log10(ratiovirushost)))+
  geom_boxplot()+ 
  geom_point(aes(colour = cast_type), size = 7,
             position = position_jitterdodge())+  
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(lineheight = 2),  
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line.x.bottom = element_line(colour = 'black', size = 0.7),
        axis.line.y.left = element_line(colour = 'black', size = 0.7),
        axis.ticks.length = unit(.18, "cm"),
        axis.ticks = element_line(colour = 'black', size = 0.7))

ser_protease <- subtilisin_exp_updated %>% subset(day == "7") %>% 
  filter(padj < 0.05 & log2FoldChange < -1) %>% 
  ggplot(aes(x=as.factor(Gene), y=log2FoldChange))+ 
  geom_bar(stat = "identity", fill = "steelblue")+
  geom_text(aes(label = log2FoldChange), vjust = 0.005, color = "white", size = 3.5)+
  theme_classic()+
  theme(legend.position = "none",
        legend.text = element_text(lineheight = 2),  
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        axis.line.x.bottom = element_line(colour = 'black', size = 0.7),
        axis.line.y.left = element_line(colour = 'black', size = 0.7),
        axis.ticks.length = unit(.18, "cm"),
        axis.ticks = element_line(colour = 'black', size = 0.7))+
  scale_y_continuous(name = "Fold change (log2)", breaks = c(0,-1,-2,-3,-4,-5))


ggarrange(virus_host, ser_protease, align = 'hv', ncol = 2, nrow = 2, legend = "bottom", common.legend = TRUE)

#serine protease and caspase activity during infection of calcified cells#########################################################################################
sub_casp <- induced_media %>% subset(experiment.name == "fluorescent_inhibitor_infection_2")

sub_casp_cells <- sub_casp %>%  summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"))
sub_casp_liths <- sub_casp %>%  summarySE(measurevar = "lith.concentration", groupvars = c("hour","treatment"))
sub_casp_virus <- sub_casp %>%  summarySE(measurevar = "virus.concentration", groupvars = c("hour","treatment"))


sub_casp$flisp_per_cell <- sub_casp$flisp_rfu / sub_casp$cell.count
sub_casp_flisp <- sub_casp %>%  summarySE(measurevar = "flisp_per_cell", groupvars = c("hour","treatment"))
sub_casp_flisp_rfu <- sub_casp %>%  summarySE(measurevar = "flisp_rfu", groupvars = c("hour","treatment"))
sub_casp_perc_pos_sum <- sub_casp %>%  summarySE(measurevar = "flisp_perc_pos", groupvars = c("hour","treatment"))

sub_casp_casp <- sub_casp %>%  summarySE(measurevar = "caspase_perc_pos", groupvars = c("hour","treatment")) 


cells <- sub_casp_cells %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "infected"))) %>%  
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment, shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120))+
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("Control", "Cells + EhVs"),
                      values = c("control" = "#000000",
                                 "infected" = "#000000"))+
  scale_shape_manual(values = c(19,21))

#source data
# write.csv(sub_casp_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs3b.csv")

liths <- sub_casp_liths %>%
  arrange(lith.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "infected"))) %>%  
  ggplot(aes(x=hour, y=log10(lith.concentration), colour=treatment, shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(lith.concentration-se), ymax=log10(lith.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120))+
scale_y_continuous(name = bquote(atop("Coccolith concentration; log10", "coccoliths mL"^-1)))+
scale_colour_manual("Treatment",
                    labels = c("Control", "Cells + EhVs"),
                    values = c("control" = "#000000",
                               "infected" = "#000000"))+
scale_shape_manual(values = c(19,21))

#source data
# write.csv(sub_casp_liths, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs3c.csv")

virus  <- sub_casp_virus %>%
  arrange(virus.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "infected"))) %>%  
  ggplot(aes(x=hour, y=log10(virus.concentration), colour=treatment, shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(virus.concentration-se), ymax=log10(virus.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120))+
  scale_y_continuous(name = bquote(atop("Virus concentration; log10", "virus mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("Control", "Cells + EhVs"),
                      values = c("control" = "#000000",
                                 "infected" = "#000000"))+
  scale_shape_manual(values = c(19,21))

#source data
# write.csv(sub_casp_virus, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs3d.csv")

sub_casp %>%
  filter(hour < 48) %>%
  arrange(flisp_perc_pos) %>%  
  mutate(treatment = factor(treatment, levels=c("control",  
                                                "infected"))) %>%  
  ggplot(aes(x=as.factor(hour), y=flisp_perc_pos, colour=treatment, shape = treatment))+  
  geom_boxplot(position = )+
  geom_point(aes(colour = treatment), size = 5,
             position = position_jitterdodge())+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.8,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
# scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+  
# scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
scale_colour_manual("Treatment",
                    labels = c("Control", "Cells + EhVs"),
                    values = c("control" = "#000000",
                               "infected" = "#000000"))+
scale_shape_manual(values = c(19,21,19,21))

#prep file
sub_casp <- sub_casp %>%
  filter(hour < 48)

#source data
# write.csv(sub_casp, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig1d.csv")

#subsetting the data into two time points
t0_flisp <- subset(sub_casp, hour == "0")
bartlett.test(flisp_perc_pos ~ treatment, data = t0_flisp)
t.test(flisp_perc_pos ~ treatment, var.equal = TRUE, conf.level = 0.95, data = t0_flisp)
#pvalue is non-significant 0.15
t24_flisp <- subset(sub_casp, hour < 25)
bartlett.test(flisp_perc_pos ~ treatment, data = t24_flisp)
t.test(flisp_perc_pos ~ treatment, var.equal = TRUE, conf.level = 0.95, data = t24_flisp)
#pvalue is significant 0.005187

ggarrange(cells, liths, virus, align = 'hv', ncol = 2, nrow = 2)

#protease enzyme activity assay########################################################################################################################################
protease_activity_growth_data <- subset(induced_media, experiment.name == "374C_infection_2" | experiment.name == "374C_infection_2_filtrate")

cell_plot <- protease_activity_growth_data %>% filter(experiment.name == "374C_infection_2") %>%
  filter(hour < 168) %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=as.factor(hour), shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Hour",
                      values = c("0" = "#88CCEE",
                                 "24" = "#CC6677",
                                 "48" = "#117733",
                                 "72" = "#332288",
                                 "96" = "#44AA99",
                                 "120" = "#661100"))+
  scale_shape_manual(labels = c("Control", "Cells + EhVs"), values = c(19,17))

#prep data
protease_activity_growth_data <- protease_activity_growth_data %>% filter(experiment.name == "374C_infection_2") %>%
  filter(hour < 168)

#source data
# write.csv(protease_activity_growth_data, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig3c.csv")

virus_plot <- protease_activity_growth_data %>% filter(experiment.name == "374C_infection_2") %>%
  filter(hour < 168) %>%
  ggplot(aes(x=hour, y=log10(virus.concentration), colour=as.factor(hour), shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+
  scale_y_continuous(name = bquote(atop("Virus concentration; log10", "viruses mL"^-1)))+
  scale_colour_manual("Hour",
                      values = c("0" = "#88CCEE",
                                 "24" = "#CC6677",
                                 "48" = "#117733",
                                 "72" = "#332288",
                                 "96" = "#44AA99",
                                 "120" = "#661100"))+
  scale_shape_manual(labels = c("Control", "Cells + EhVs"), values = c(19,17))


bacteria_plot <- protease_activity_growth_data %>% filter(experiment.name == "374C_infection_2") %>%
  filter(hour < 168) %>%
  ggplot(aes(x=hour, y=log10(bacterial.concentration), colour=as.factor(hour), shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+
  scale_y_continuous(name = bquote(atop("Bacterial concentration; log10", "bacteria mL"^-1)))+
  scale_colour_manual("Hour",
                      values = c("0" = "#88CCEE",
                                 "24" = "#CC6677",
                                 "48" = "#117733",
                                 "72" = "#332288",
                                 "96" = "#44AA99",
                                 "120" = "#661100"))+
  scale_shape_manual(labels = c("Control", "Cells + EhVs"), values = c(19,17))

calcified_plot <- protease_activity_growth_data %>% filter(experiment.name == "374C_infection_2") %>%
  filter(hour < 168) %>%
  ggplot(aes(x=hour, y=naked, colour=as.factor(hour), shape = treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+
  scale_y_continuous(name = bquote(atop("Percentage of naked cells")))+
  scale_colour_manual("Hour",
                      values = c("0" = "#88CCEE",
                                 "24" = "#CC6677",
                                 "48" = "#117733",
                                 "72" = "#332288",
                                 "96" = "#44AA99",
                                 "120" = "#661100"))+
  scale_shape_manual(labels = c("Control", "Cells + EhVs"), values = c(19,17))

extracellular_protease <- protease_activity %>% 
  filter(hour < 168) %>%
  ggplot(aes(x=as.factor(hour), y=log10(hydrolysis_rate_nmol_min_ml), colour=as.factor(hour), shape = treatment))+
  geom_boxplot()+
  geom_point(aes(colour = as.factor(hour)),  
             size = 3, 
             position = position_jitterdodge())+
  theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+
  scale_y_continuous(name = bquote(atop("Extracellular protrease activity; log10", (nmol ~protein~min^-1~ml^-1))))+
  scale_colour_manual("Hour",
                      values = c("0" = "#88CCEE",
                                 "24" = "#CC6677",
                                 "48" = "#117733",
                                 "72" = "#332288",
                                 "96" = "#44AA99",
                                 "120" = "#661100"))+
  scale_shape_manual(labels = c("Control", "Cells + EhVs"), values = c(19,17))

#96 h
enzyme_activity_96h <- subset(protease_activity, hour == "96")
bartlett.test(hydrolysis_rate_nmol_min_ml ~ treatment, data=enzyme_activity_96h)
t.test(hydrolysis_rate_nmol_min_ml ~ treatment, data=enzyme_activity_96h)
#120 h
enzyme_activity_120h <- subset(protease_activity, hour == "120")
bartlett.test(hydrolysis_rate_nmol_min_ml ~ treatment, data=enzyme_activity_120h)
t.test(hydrolysis_rate_nmol_min_ml ~ treatment, data=enzyme_activity_120h, var.equal=FALSE, conf.level=0.95)

#cells growing in the 0.02 filtrates that I did protease activity on
protease_activity_growth_data_sum <- protease_activity_growth_data %>%  filter(experiment.name == "374C_infection_2_filtrate") %>% summarySE(measurevar = "cell.concentration", groupvars = c("hour","treatment"), na.rm = TRUE)

cell_filtrate_plot <- protease_activity_growth_data_sum %>%
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment))+  
  geom_point(size=5)+theme_bw()+  
  theme_classic()+
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Hour", breaks=c(0,24,48,72,96,120,144,168))+
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
scale_colour_manual("Hour",
                    values = c("control" = "#000000",
                               "t0_filtrate" = "#88CCEE",
                               "24h_filtrate" = "#CC6677",
                               "48h_filtrate" = "#117733",
                               "72h_filtrate" = "#332288",
                               "96h_filtrate" = "#44AA99",
                               "120h_filtrate" = "#661100"))

#source data
# write.csv(protease_activity_growth_data_sum, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/fig3b.csv")

enzyme_growth_data_96h <- subset(protease_activity_growth_data, hour == "96")
enzyme_growth_data_96h <- enzyme_growth_data_96h %>% filter(filtration_method == "anotop")
enzyme_growth_data_96h = mutate(enzyme_growth_data_96h, treatment = factor(treatment, levels=unique(treatment)))

#log10 in the model to remove heteroscedasticity
model = lm(log10(cell.concentration) ~ treatment, data = enzyme_growth_data_96h)
histogram(residuals(model), col="darkgrey")#should have normal distribution
plot(fitted(model), residuals(model))
Anova(model, type = "II") #p-value: 1.666e-07
(HSD.test(model, "treatment"))

ggarrange(extracellular_protease, cell_filtrate_plot, nrow = 2, ncol = 2, align = "hv", common.legend = TRUE, legend =  "top")

ggarrange(cell_plot, calcified_plot, virus_plot, bacteria_plot,  nrow = 2, ncol = 2, align = "hv", legend = "none")

#Does removing the bacteria delay infection?#######################################################################################################################################
#fold change in bacterial reduction as the cell are treated with antibiotcs
bacteria_fold_change <- antibiotic_bacterial_concentrations %>% 
  ggplot(aes(x = time_point, y = fold_reduction))+
  geom_bar(stat = "identity", width = 0.5)+
  theme_bw()+
  theme_classic()+ 
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  scale_y_reverse(name = "Fold reduction in bacteria")+
  scale_x_continuous(name = "Sampling time point")
  

#cell concentration
antibiotic_infection_data <- induced_media %>% filter(experiment.name == "374C_antibiotic_infection")

antibiotic_infection_cells <- antibiotic_infection_data %>% 
  summarySE(measurevar = "cell.concentration",  
            groupvars = c("hour","strain_id","treatment"))

antibiotic_infection_data_96h <- antibiotic_infection_data %>% filter(hour > 95)

#cell concentrations for the antibiotic infection as a boxplot
boxplot <- antibiotic_infection_data_96h %>% 
  subset(treatment == "infected" |treatment == "infected_antibiotics") %>%
  arrange(cell.concentration) %>%
  mutate(treatment = factor(treatment, levels=c("control_heat_killed",
                                                "infected",
                                                "control_antibiotics_heat_killed",
                                                "infected_antibiotics"))) %>%
  ggplot(aes(x=as.factor(treatment), y=cell.concentration, fill = treatment)) + 
  geom_boxplot(alpha=0.2)+
  geom_point(aes(colour = treatment), size = 5,
             position = position_jitterdodge())+  
  theme_classic()+ 
  theme(legend.position = "none",  
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),  
        title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  labs(x = "", y ="")+
  scale_y_continuous(limits = c(900000,1350000))+
  scale_fill_manual("Treatment",
                       labels = c("Cells + heat-killed EhVs","Cells + EhVs","Cells + heat-killed EhVs", "Cells + EhVs"),
                       values = c("control_heat_killed" = "#FFCC00",
                                  "infected" = "#FFCC00",
                                  "control_antibiotics_heat_killed" = "#97cfd0",
                                  "infected_antibiotics" = "#97cfd0"))+
  scale_color_manual("Treatment",
                    labels = c("Cells + heat-killed EhVs","Cells + EhVs","Cells + heat-killed EhVs", "Cells + EhVs"),
                    values = c("control_heat_killed" = "#FFCC00",
                               "infected" = "#FFCC00",
                               "control_antibiotics_heat_killed" = "#97cfd0",
                               "infected_antibiotics" = "#97cfd0"))

#plotting cell concentrations for the antibiotic infection
cells <- antibiotic_infection_cells %>%
  subset(treatment == "infected" |treatment == "infected_antibiotics") %>%
  # arrange(cell.concentration) %>%
  # mutate(treatment = factor(treatment, levels=c("control",  
  #                                               "viral im",
  #                                               "30kDA", 
  #                                               "10kDA",
  #                                               "3kDA"))) %>%  
  ggplot(aes(x=hour, y=log10(cell.concentration), colour=treatment))+
  geom_point(size=5)+theme_bw()+
  theme_classic()+ 
  theme(legend.position = "none",  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=log10(cell.concentration-se), ymax=log10(cell.concentration+se)), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192))+  
  scale_y_continuous(name = bquote(atop("Cell concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("Cells + heat-killed EhVs","Cells + EhVs","Cells + heat-killed EhVs", "Cells + EhVs"),
                      values = c("control_heat_killed" = "#FFCC00",
                                 "infected" = "#FFCC00",
                                 "control_antibiotics_heat_killed" = "#97cfd0",
                                 "infected_antibiotics" = "#97cfd0"))

#source data
# write.csv(antibiotic_infection_cells, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs8b.csv")

antibiotic_cells <- cells + inset_element(boxplot, left = 0.4, bottom = 0.001, right = 1, top = 0.55)
antibiotic_cells <- as.ggplot(antibiotic_cells)

#t.test for the cell concentrations at 96 h between the infected treatments with and without antibiotics
antibiotic_infection_data_96h_infected <- antibiotic_infection_data_96h %>% 
  subset(treatment == "infected" |treatment == "infected_antibiotics")
t.test(log10(cell.concentration) ~ treatment, data = antibiotic_infection_data_96h_infected)

#virus concentration
antibiotic_virus <- antibiotic_infection_data %>%
  subset(treatment == "infected" |treatment == "infected_antibiotics") %>%
  ggplot(aes(x=hour, y=log10(virus.concentration), colour = treatment))+  
  geom_line(data=antibiotic_infection_data[!is.na(antibiotic_infection_data$virus.concentration),],aes(colour=treatment,group=treatment))+
  geom_point(size = 5)+
  theme_bw()+
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  scale_x_continuous(name = "Time (h)", limits = c(0,96), breaks=c(0,24,48,72,96,120,144,168,192))+  
  scale_y_continuous(name = bquote(atop("Virus concentration; log10", "cells mL"^-1)))+
  scale_colour_manual("Treatment",
                      labels = c("Cells + heat-killed EhVs","Cells + EhVs","Cells + heat-killed EhVs", "Cells + EhVs"),
                      values = c("control_heat_killed" = "#FFCC00",
                                 "infected" = "#FFCC00",
                                 "control_antibiotics_heat_killed" = "#97cfd0",
                                 "infected_antibiotics" = "#97cfd0"))

#source data
# write.csv(antibiotic_infection_data, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs8d.csv")

#percentage of naked cells
antibiotic_infection_naked <- antibiotic_infection_data %>% summarySE(measurevar = "naked",  
            groupvars = c("hour","treatment"))

antibiotic_naked <- antibiotic_infection_naked %>%
  subset(treatment == "infected" |treatment == "infected_antibiotics") %>%
  # arrange(cell.concentration) %>%
  # mutate(treatment = factor(treatment, levels=c("control",  
  #                                               "viral im",
  #                                               "30kDA", 
  #                                               "10kDA",
  #                                               "3kDA"))) %>%  
  ggplot(aes(x=hour, y=naked, colour=treatment))+  
  geom_point(size=5)+theme_bw()+ 
  theme_classic()+  
  theme(legend.position = c(0.1,0.88), legend.text = element_text(lineheight = 2),  
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),  
        axis.title=element_text(size=16),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        axis.line = element_line(colour = 'black', size = 0.75),
        axis.ticks.length=unit(.15, "cm"),  
        axis.ticks = element_line(colour = "black", size = 1))+  
  geom_errorbar(aes(ymin=naked-se, ymax=naked+se), width= 5, size=0.75)+  
  geom_line(aes(colour=treatment,group=treatment))+
  scale_x_continuous(name = "Time (h)", breaks=c(0,24,48,72,96,120,144,168,192))+  
  scale_y_continuous(name = "Percentage of naked cells")+
  scale_colour_manual("Treatment",
                      labels = c("Cells + heat-killed EhVs","Cells + EhVs","Cells + heat-killed EhVs", "Cells + EhVs"),
                      values = c("control_heat_killed" = "#FFCC00",
                                 "infected" = "#FFCC00",
                                 "control_antibiotics_heat_killed" = "#97cfd0",
                                 "infected_antibiotics" = "#97cfd0"))

#source data
# write.csv(antibiotic_infection_naked, "/Users/Christopher/Dropbox/Lith Adsorption/Spreadsheets/induced_media/exported_source_data_files/figs8c.csv")

#residuals for the comparison of naked cells looked ok
t.test(naked ~ treatment, data = antibiotic_infection_data_96h_infected)

#plotting everything together
ggarrange(bacteria_fold_change, antibiotic_cells, antibiotic_naked, antibiotic_virus, nrow = 2, ncol = 2, common.legend = TRUE)
