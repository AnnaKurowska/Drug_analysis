##Packages
library(tidyr)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(tsne)
library(gridExtra)

#useful functions
#everything()
#any

##Loading the data
setwd("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical")
data_trial <- read.csv("Drugs_trial.csv")

## Conditions, which drugs induce desired gene expression
#modifying data frame

neutral_genes <-"BM2"
anti_genes <- c("HBD1",	"HBD2",	"HBD3",	"HBD4",	"LL37",	"REG3A",	"PLA2G2A")
pro_genes <- c("IL1B","IL8","TNF")

anti_sums <- c()  ## creating lists inside
pro_sums <- c() 
new_data <- list()
for (i in 1:nrow(data_trial)) {
  if ((max(data_trial[i,anti_genes]) >= 2) & (max(data_trial[i,pro_genes]) < 10)){
    new_data[[i]] <- data_trial[i,]
    anti_sums[i] <- (sum(data_trial[i,anti_genes]))/(length(data_trial[i,anti_genes]))
    pro_sums[i] <- (sum(data_trial[i,pro_genes]))/(length(data_trial[i, pro_genes]))
  }
}

df <-as.data.frame(do.call("rbind",new_data))
ratios <-anti_sums/pro_sums
ratios <- as.data.frame(ratios[!is.na(ratios)])
ratios$Drug <- df$Drug
names(ratios) <- c("Ratio", "Drug")
df <-full_join(df,ratios, by = "Drug")
##### maybe make that into one function!


# representing it as barplots for each compound separately
df_plotting <-gather(df, key = "Gene", value = "ExpressionRatio", -c(Drug,Ratio))

potential_drugs <- df_plotting$Drug

estradiol <-df_plotting[df_plotting$Drug == "Estradiol",]
estradiol$Gene <- factor(estradiol$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNF"))

ggplot(estradiol, aes(x = Gene, y = ExpressionRatio)) +
  geom_bar(stat = "identity") +
  theme_bw()
#...
#### Representing the ratios
ggplot(df_plotting, aes(x = Drug, y = Ratio, fill = Drug)) + 
  geom_bar(stat = "identity") + 
  theme_classic()
 ## the bigger the ratio, the better: means that higher average expression of 
 ##could also just do the sum, bc ratio may be wrong if only 1 inflammatory gene is on...
### Version with the sum instead of average expressions (+ as a function): 

represent_data <- function(dataset){
  
sum_ratio <- c()  ## creating lists inside <- c() 
new_data <- list()

for (i in 1:nrow(dataset)) {
 if ((max(dataset[i,anti_genes]) >= 2) & (max(dataset[i,pro_genes]) < 10)) {
  new_data[[i]] <- dataset[i,]
  sum_ratio[i] <- (sum(dataset[i,anti_genes]))/ (sum(dataset[i,pro_genes]))
 }
}
  df <-as.data.frame(do.call("rbind",new_data))
  ratios <- as.data.frame(sum_ratio[!is.na(sum_ratio)])
  ratios$Drug <- df$Drug
  names(ratios) <- c("Ratio", "Drug")
  df2 <-full_join(df,ratios, by = "Drug")
  return(df2)
}

#call function
see <-represent_data(data_trial)
see_plotting <-gather(see, key = "Gene", value = "ExpressionRatio", -c(Drug,Ratio))

#plot
ggplot(see_plotting, aes(x = Drug, y = Ratio, fill = Drug)) + 
  geom_bar(stat = "identity") + 
  theme_classic()

### now data with different concentrations 
data_conc <- read.csv("drugs_concentration.csv")

## firs replace dilutions into 'No_dil' and 'dil' until you know real concentrations

# d: digit, s:space w:character
for (i in 1:nrow(data_conc)){ 
  
  if (str_detect(data_conc[i,1], "\\d\\s\\w+") == TRUE){ 
    data_conc[i,1] <- "Dil"
    
  } else {
    data_conc[i,1] <- "No_Dil"
  }
}

### choose both rows for drugs that match the conditions at least once
 #tutaj bez znajdowania zadnych values, tylko for visualization
  c_data <- list()
  for (i in 1:nrow(data_conc)) {
    
    if ((max(data_conc[i,anti_genes]) >= 2) & (max(data_conc[i,pro_genes]) < 10)) {
      if (str_detect(data_conc[i,1], "No") == TRUE){
        
        c_data[[i]] <- data_conc[i,]
        c_data[[i+1]] <- data_conc[i+1,]
        
      } else{
        
        c_data[[i]] <- data_conc[i,]
        c_data[[i-1]] <- data_conc[i-1,]
      }
    }
  }

##Now visualize the data
df_c <-as.data.frame(do.call("rbind",c_data)) # data after selection of drugs that met both conditions
c_modified <- gather(df_c, key = "Gene", value = "Expression", -c(Condition, Drug))

########################## Genevestigator data ############################

# Main points about the data:
# 1. Different amount of conditions in each case, but treated as a separate sample
# 2. Sometimes multiple: cell lines, concentrations, exposure times.

###idea is, find compound names that satisfy the conditions and then include them in a new table by selecting them 

##Identifying drugs that meet conditions (a way more efficient approach)
genes_names <- c()
for (i in 1:nrow(data_conc)) {
  
  if ((max(data_conc[i,anti_genes]) >= 2) & (max(data_conc[i,pro_genes]) < 10)) {
   genes_names[i] <- data_conc[i,2]
   genes_names <- genes_names[!is.na(genes_names)]
   genes_names <- genes_names[!duplicated(genes_names)]
  }
}

##Now only select rows that match your compounds that met the conditions
concentration_data <-data_conc[data_conc$Drug %in% genes_names,]
   # and now you can continue with the analysis, depending on what they do

### EXERCISES with the full table that i copied manually###################
setwd("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical")
full <- read.csv("full_table.csv")

neutral_genes <-"BM2"
anti_genes <- c("HBD1",	"HBD2",	"HBD3",	"HBD4",	"LL37",	"REG3A",	"PLA2G2A")
pro_genes <- c("IL1B","IL8","TNFa")

for (i in 4:nrow(full)){ 
  
  if (str_detect(full[i,1], "\\d\\s\\w+") == TRUE){ 
    full[i,1] <- "10 um"
    
  } else {
    full[i,1] <- "100 um"
  }
}

compound_names <- c()
for (i in 4:nrow(full)) {
  
  if ((max(full[i,anti_genes]) >= 2) & (max(full[i,pro_genes]) < 10)) {
    compound_names[i] <- full[i,2]
    compound_names <- compound_names[!is.na(compound_names)]
    compound_names <- compound_names[!duplicated(compound_names)]
  }
}

concentration_data <-full[full$Nom.Echantillon %in% compound_names,] ###  filtered data according to conditions
## representing data
length(compound_names) ### too many of them to represent them all


### Identify compounds with aberrant B2M ratio!##########
concentration_data
compound_names_B2M <- c()
for (i in 1:nrow(concentration_data)) {
  
  if ((max(concentration_data[i, "B2M"]) < 0.5) | (max(concentration_data[i,"B2M"]) > 2)) {
    compound_names_B2M[i] <- concentration_data[i,2]
    compound_names_B2M <- compound_names_B2M[!is.na(compound_names_B2M)]
    compound_names_B2M <- compound_names_B2M[!duplicated(compound_names_B2M)]
  }
}
# so flagged compounds!!:
compound_names_B2M

#calculate compounds with the greatest ratios
ratios <- c()
for (i in 1:nrow(concentration_data)){
  ratios[i] <- (sum(concentration_data[i,anti_genes]))/ (sum(concentration_data[i,pro_genes]))
}
concentration_data$Ratio <-ratios
  #represent compounds in the order of the greatest
ordered_C_data <-concentration_data[with(concentration_data,order(-Ratio)),]
write.csv(head(ordered_C_data), "head_ratio_values.csv")

  #VISUALISE RATIOS (remember NO COMPOUND LABELS)########### change labels to micro!

ratio_data_full_table <-concentration_data %>%
  select(Conditions,Nom.Echantillon, Ratio) 

ggplot(ratio_data_full_table, aes(x = Nom.Echantillon, y = Ratio, fill = Conditions)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  coord_flip() +
  theme(axis.text.y=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  xlab("Compounds")
  



##here, regardless of the ratio, do all compounds for each gene expression graph. As suggested by Brice (a scaterplot)

#Hits: represent all data! From "full" dataset
no_dil_data <- full[full$Conditions == "No_Dil",]

#HBD1
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = HBD1, colour = dplyr::case_when(HBD1 >= 2 ~ ">2", HBD1 < 2 ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("HBD1 Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 4))
  
HBD1_filtered_genes <- no_dil_data[no_dil_data$HBD1 >= 2,2]

#HBD2
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = HBD2, colour = dplyr::case_when(HBD2 >= 2 ~ ">2", HBD2 < 2 ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("HBD2 Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(trans = "log2")

HBD2_filtered_genes <- no_dil_data[no_dil_data$HBD2 >= 2,2]

#HBD3
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = HBD3, colour = dplyr::case_when(HBD3 >= 2 ~ ">2", HBD3 < 2 ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("HBD3 Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(trans = "log2")

HBD3_filtered_genes <- no_dil_data[no_dil_data$HBD3 >= 2,2]

intersect(intersect(HBD1_filtered_genes, HBD2_filtered_genes), HBD3_filtered_genes)


#HBD4
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = log2(HBD4), colour = dplyr::case_when(log2(HBD4) >= log2(2) ~ ">2", log2(HBD4) < log2(2) ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("HBD4 Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 9))

no_dil_data[no_dil_data$HBD4 >= 2,2]

#LL37
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = log2(LL37), colour = dplyr::case_when(log2(LL37) >= log2(2) ~ ">2", log2(LL37) < log2(2) ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("LL37 Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(-2, 6))

no_dil_data[no_dil_data$LL37 >= 2,2]

#REG3A
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = log2(REG3A), colour = dplyr::case_when(log2(REG3A) >= log2(2) ~ ">2", log2(REG3A) < log2(2) ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("REG3A Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
  
no_dil_data[no_dil_data$REG3A >= 2,2]

#PLA2G2A
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = log2(PLA2G2A), colour = dplyr::case_when(log2(PLA2G2A) >= log2(2) ~ ">2", log2(PLA2G2A) < log2(2) ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("PLA2G2A Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

no_dil_data[no_dil_data$PLA2G2A >= 2,2]

### Check the same for pro-inflammatory genes

#IL1B
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = IL1B, colour = dplyr::case_when(IL1B >= 2 ~ ">2", IL1B < 2 ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("IL1B Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
  
#IL8
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = log2(IL8), colour = dplyr::case_when(log2(IL8) >= log2(2) ~ ">2", log2(IL8) < log2(2) ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("IL8 Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

#TNFa
ggplot(no_dil_data, aes(x = Nom.Echantillon, y = log2(TNFa), colour = dplyr::case_when(log2(TNFa) >= log2(2) ~ ">2", log2(TNFa) < log2(2) ~ "<2 "))) +
  geom_point() +
  theme_bw() + 
  scale_colour_discrete("Activity") +
  xlab("Compounds") +
  ylab("TNFa Expression Ratio") +
  guides(size = "none") +
  # geom_segment(aes(x=Nom.Echantillon, xend=Nom.Echantillon, y=0, yend=HBD1)) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
  
### Representation of ordered compounds 
compounds <-unique(full[4:nrow(full),"Nom.Echantillon"])
ordered_compounds <-c("2-MethoxyestriadIol","Vincistrine sulfate","Celastrol","Progesterone","Estriadol","Oxytetracyline (Terramycin)","4-Demethylepipodophyllotoxin ","Colchicine","Nocodazole","Dehydrocostus Lactone","Kitasamycin")

gathered_data<- gather(full[4:nrow(full),], key = "Gene", value = "Expression", -c(Conditions,Nom.Echantillon))

#MethoxyestriadIol
MethoxyestriadIol <-gathered_data[gathered_data$Nom.Echantillon == "2-MethoxyestriadIol",]
MethoxyestriadIol$Gene <- factor(MethoxyestriadIol$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
MethoxyestriadIol$Conditions <- factor(MethoxyestriadIol$Conditions, levels = c("No_Dil", "Dil"))

ggplot(MethoxyestriadIol, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) 
 
#Vincistrine sulfate

Vincistrine <-gathered_data[gathered_data$Nom.Echantillon == "Vincistrine sulfate",]
Vincistrine$Gene <- factor(Vincistrine$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Vincistrine$Conditions <- factor(Vincistrine$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Vincistrine, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM"))

#Celastrol
Celastrol <-gathered_data[gathered_data$Nom.Echantillon == "Celastrol",]
Celastrol$Gene <- factor(Celastrol$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Celastrol$Conditions <- factor(Celastrol$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Celastrol, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, NA)) 

#Progesterone
Progesterone <-gathered_data[gathered_data$Nom.Echantillon == "Progesterone",]
Progesterone$Gene <- factor(Progesterone$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Progesterone$Conditions <- factor(Progesterone$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Progesterone, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) 

#Estriadol
Estriadol <-gathered_data[gathered_data$Nom.Echantillon == "Estriadol",]
Estriadol$Gene <- factor(Estriadol$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Estriadol$Conditions <- factor(Estriadol$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Estriadol, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) 

#Oxytetracyline (Terramycin)

Oxytetracyline <-gathered_data[gathered_data$Nom.Echantillon == "Oxytetracyline (Terramycin)",]
Oxytetracyline$Gene <- factor(Oxytetracyline$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Oxytetracyline$Conditions <- factor(Oxytetracyline$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Oxytetracyline, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, NA)) 

#4-Demethylepipodophyllotoxin

Demethylepipodophyllotoxin <-gathered_data[gathered_data$Nom.Echantillon == "4-Demethylepipodophyllotoxin ",]
Demethylepipodophyllotoxin$Gene <- factor(Demethylepipodophyllotoxin$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Demethylepipodophyllotoxin$Conditions <- factor(Demethylepipodophyllotoxin$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Demethylepipodophyllotoxin, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) +
  scale_y_continuous(expand = c(0, 0),limits = c(0, NA)) 

#Colchicine

Colchicine <-gathered_data[gathered_data$Nom.Echantillon == "Colchicine",]
Colchicine$Gene <- factor(Colchicine$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Colchicine$Conditions <- factor(Colchicine$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Colchicine, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM"))

#Nocodazole
Nocodazole <-gathered_data[gathered_data$Nom.Echantillon == "Nocodazole",]
#filter(flights, month == 1, day == 1)
Nocodazole$Gene <- factor(Nocodazole$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Nocodazole$Conditions <- factor(Nocodazole$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Nocodazole, aes(x = Gene, y = log2(Expression), fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) 

#Dehydrocostus Lactone

Dehydrocostus <-gathered_data[gathered_data$Nom.Echantillon == "Dehydrocostus Lactone",]
Dehydrocostus$Gene <- factor(Dehydrocostus$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Dehydrocostus$Conditions <- factor(Dehydrocostus$Conditions, levels = c("No_Dil", "Dil"))

ggplot(Dehydrocostus, aes(x = Gene, y = Expression, fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) +
  scale_y_continuous(trans='log2')

#Kitasamycin
Kitasamycin <-gathered_data[gathered_data$Nom.Echantillon == "Kitasamycin",]
#filter(flights, month == 1, day == 1)
Kitasamycin$Gene <- factor(Kitasamycin$Gene, levels = c("B2M","HBD1","HBD2", "HBD3","HBD4","LL37","REG3A","PLA2G2A","IL1B","IL8","TNFa"))
Kitasamycin$Conditions <- factor(Kitasamycin$Conditions, levels = c("No_Dil", "Dil"))
Kitasamycin$Start <- 0.1

ggplot(Kitasamycin, aes(x = Gene, y = Expression, fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM"))  +
  scale_y_continuous(trans='log2') ## issue bc you have no log of 0
  


ggplot(Kitasamycin, aes(x = Gene, y = Expression, fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  ylab(label = "Expression Ratio") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_discrete(name = "Concentration", labels=c("100 μM","10 μM")) +
  scale_y_continuous(trans='log2') ## issue bc you have no log of 0 +

write.csv(Kitasamycin,"Kitasamycin.csv")

#alternatively a table for compounds with the greatest ratio
most_DE_compounds <- c("Colchicine", "4-Demethylepipodophyllotoxin", "Genistin", "Dehydrocostus Lactone", "Oxytetracyline (Terramycin)","Vincistrine sulfate")

most_DE_table <-concentration_data[concentration_data$Nom.Echantillon %in% most_DE_compounds,]
write.csv(most_DE_table, "most_DE_compounds.csv") # save it as an output for representation

## visualize expression ratio for all drugs
gathered_data
ggplot(gathered_data, aes(x = Nom.Echantillon, y = Ratio, fill = Conditions)) +
  geom_bar(stat = "Identity",position = "dodge") + 
  coord_flip() +
  theme_classic() +
  xlab("Compound") +
  ylab(" Ratio (DE AMP/ DE Inflammatory)") 
  ## remember these are data that follow your conditions


#compute statistics for each gene across studies:
apply(concentration_data[,3:13], 2, mean)
apply(concentration_data[,3:13], 2, median)
apply(concentration_data[,3:13], 2, max)
apply(concentration_data[,3:13], 2, min)

### Since 37 conditions match my initial search, maybe i need to do do cross-validation

# let's start with compounds from the most most_DE_table
most_DE_compounds

 ###################CROSS VALIDATION OF MY TEST COMPOUNDS ###########
#1. I'll start with Genevestigator and if not GEO

## data for ISOFLAVONE (genistin) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60971. Just remember its not actually genistin but genistein

try_df <- read.delim("GSE60971_normalized.txt")

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
library(BiocManager)
BiocManager::install("illuminaHumanv4.db")
library(illuminaHumanv4.db)

x <- illuminaHumanv4ALIAS2PROBE 
mapped_probes <- mappedkeys(x)
xx <- as.list(x[mapped_probes])
genes_id <- c("BM2" ,    "HBD1" ,   "HBD2" ,   "HBD3"  ,  "HBD4" ,   "LL37",  "REG3A" , "PLA2G2A", "IL1B"  ,  "IL8"    , "TNF" )
genes_with_IDs <-do.call("rbind",xx[genes_id])
genes_to_retrieve <- c("ILMN_1686573","ILMN_2132515","ILMN_1688580","ILMN_1757504","ILMN_2382679","ILMN_1740586","ILMN_1775501","ILMN_1666733", "ILMN_2184373","ILMN_1728106")

genistin_crossdata <- try_df[try_df$ID_REF %in% genes_to_retrieve,] 
# genistin_crossdata <-as.data.frame(t(genistin_crossdata))
gene_names <- c("IL8", "HBD1", "LL37", "TNF", "PLA2G2A","REG3A", "IL1B", "HBD3", "IL8_2", "REG3A_2")
new_colnames_genistin <- c("ID_REF", "24h30_1", "24h60_1", "24h100_1", "24hK_1", "48h30_1", "48h60_1", "48h100_1", "48hK_1", "24h30_2", "24h60_2", "24h100_2", "24hK_2", "48h30_2", "48h60_2", "48h100_2", "48hK_2", "24h30_3", "24h60_3", "24h100_3", "24hK_3", "48h30_3", "48h60_3", "48h100_3", "48hK_3")
colnames(genistin_crossdata) <- new_colnames_genistin
genistin_crossdata <- genistin_crossdata[,1:9]
genistin_crossdata$Genes  <- gene_names
gathered_genistin <-gather(genistin_crossdata, key = "Condition",value = DE, -c(ID_REF,Genes))
gathered_genistin$Condition <- factor(gathered_genistin$Condition, levels = new_colnames_genistin )

ggplot(gathered_genistin, aes(x= Genes, y = Condition, fill = DE)) +
  geom_tile(stat = "identity") + 
  theme_classic()
### no takie mieszane wyniki bo z jednej strony HBD1 wzrasta ale jest przygłuszone pro-inflammatory genes
  # pracują nad komórkami z chorobą, dlatego widzimy taki high pro-inflammatory level


## Colchicine https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE4090
colchicine_table <- read.delim("GSE4090_series_matrix.txt")
#change column
new_colnames_colchicine_table <- c("ID_REF", "HUVEC exposure no colchicine 0_2", "HUVEC exposure 1ug/ml colchicine 24 hr", "HUVEC exposure 1 ug/ml colchicine 12 hr_1", "HUVEC exposure no colchicine 0_3", "HUVEC exposure 1 ug/ml colchicine 120'_2", "HUVEC exposure 1 ug/ml colchicine 12hr_2", "HUVEC exposure no colchicine 12 hr_1", "HUVEC exposure 1ug/ml colchicine 30'_1", "HUVEC exposure 1ug/ml colchicine 30'_2", "HUVEC exposure no colchicine 12 hr_2", "HUVEC exposure no colchicine 24 hr_1", "HUVEC exposure 100 ng/ml colchicine 120'", "HUVEC exposure 1 ug/ml colchicine 120'_1", "HUVEC exposure 100 ng/ml colchicine 12 hr", "HUVEC exposure no colchicine 0_1", "HUVEC exposure no colchicine 24 hr_2")
colnames(colchicine_table) <-new_colnames_colchicine_table

# identifying genes to probes
Vincistrine_table <- read.delim("GSE162789_series_matrix.txt")
Vincistrine_selected <- Vincistrine_table[,c(1,26,28,35,37)] # only selecting relevant samples, control + vincristine treatment
Vincistrine_names <- c("ID","EW7 + control", "EW7 + Vincristine", "SK-N-MC + control", "SK-N-MC + Vincristine")
colnames(Vincistrine_selected) <- Vincistrine_names

BiocManager::install("hugene10stprobeset.db")
library(hugene10stprobeset.db)
xx2 <- as.list(hugene10stprobesetALIAS2PROBE)
xx2[genes_id]
genes_to_retrieve <- c("8149100","8149101","8149102","8149103","8079591","8079592","8079593","8079594", "8079595","8079596","8079597","8053342","8053343", "8053344", "8053345", "8053346", "7913217", "7913218", "7913219", "7913220", "7913221", "7913222", "7913223","7913224", "7913225", "8054723", "8054724", "8054725" ,"8054726" ,"8054727" ,"8054728", "8054729", "8095681", "8095682", "8095683", "8095684", "8095685", "8095686", "8095687","8118143", "8118144", "8118145" ,"8118146" ,"8177984", "8177985" , "8177986", "8177987", "8179264", "8179265", "8179266", "8179267")

Vinc_crossdata <- Vincistrine_selected[Vincistrine_selected$ID %in% genes_to_retrieve,] 
xx2["TNF"]

########## END OF CROSS VALIDATION
