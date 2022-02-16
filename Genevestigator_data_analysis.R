library(readxl)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(hrbrthemes) # for heatmaps
library(magrittr)
library(pheatmap)
library(MixGHD) # hierarchical clustering which i will finally not use
library(viridis)
library(ggpubr)


### THe pourpose is to first filter through each of the 4 platforms to see which AMP > 2 and PRO < 10. Once I have that, I represent the ratios for only one study of each drug (sometimes they are tested on multiple cell lines/ different concentrations)

############# Automating the code
DE_Genevestigator <- function(data, AMP_data, PRO_data, AMP = 2, PRO = 10){
  
  compounds_to_be_stored <- c()
  
  for (i in 1:nrow(data)) {
    if ((max(data[i,AMP_data]) >= AMP) & (max(data[i,PRO_data]) < PRO)) {
      
      compounds_to_be_stored[i] <- data[i,1]
      
    }
  }
  compounds_to_be_stored <- compounds_to_be_stored[!is.na(compounds_to_be_stored)]
  compounds_to_be_stored <- compounds_to_be_stored[!duplicated(compounds_to_be_stored)]
  
  return(data[data$Perturbations %in% compounds_to_be_stored,])
}


###U133A_2 dataset

U133A_2 <- read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data/Toxicology-HS_AFFY_U133A-2.xlsx")

#modifying the table for analysis 
U133A_2 <-U133A_2[,-c(2:5)]
colnames(U133A_2) <- c("Perturbations", "log2:BM2", "FC:BM2","p-value:BM2",
                       "log2:DEFB1", "FC:DEFB1","p-value:DEFB1",
                       "log2:DEFB4A", "FC:DEFB4A","p-value:DEFB4A",
                       "log2:CAMP", "FC:CAMP","p-value:CAMP",
                       "log2:REG3A", "FC:REG3A","p-value:REG3A",
                       "log2:PLA2G2A", "FC:PLA2G2A","p-value:PLA2G2A",
                       "log2:IL1B", "FC:IL1B","p-value:IL1B",
                       "log2:CXCL8", "FC:CXCL8","p-value:CXCL8",
                       "log2:TNF", "FC:TNF","p-value:TNF" )
            #not all genes were found on this platform!
U133A_2 <- U133A_2[5:nrow(U133A_2),]

# turn all character columns into numeric 
class(U133A_2)
U133A_2_modified <-data.frame(lapply(U133A_2,as.numeric))
U133A_2_modified$Perturbations <- U133A_2$Perturbations
str(U133A_2_modified)

##identify compounds that meet your conditions AMP >2, PRO <20

# first define genes in AMP and PRO for which you'll run the loop
AMP_U133A_2 <- c("FC.BM2","FC.DEFB1","FC.DEFB4A", "FC.CAMP","FC.REG3A","FC.PLA2G2A")
PRO_U133A_2 <- c("FC.IL1B","FC.CXCL8","FC.TNF")

results_U133A <- DE_Genevestigator(U133A_2_modified,AMP_data = AMP_U133A_2,PRO_data = PRO_U133A_2)

finding_selected_compounds <-function(results_data){
  selected_words <- c()
  for (i in 1:length(results_data$Perturbations)){
    selected_words[i] <- str_extract(results_data$Perturbations[i], "([^\\s]+)") }
    return(unique(selected_words))
}

length(finding_selected_compounds(results_U133A))

selected_words_U133 <- c()
for (i in 1:length(results_U133A$Perturbations)){
  selected_words_U133[i] <- str_extract(results_U133A$Perturbations[i], "([^\\s]+)") 
}
selected_words_U133 <- selected_words_U133[!duplicated(selected_words_U133)]
length(selected_words_U133)

###U133PLUS_2-0 dataset

U133PLUS_2_0 <- read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data/Toxicology-HS_AFFY_U133PLUS_2-0.xlsx")

#modifying the table for analysis 
U133PLUS_2_0 <-U133PLUS_2_0[,-c(1,3,4,5)]
colnames(U133PLUS_2_0) <- c("Perturbations", "log2:BM2", "FC:BM2","p-value:BM2",
                       "log2:DEFB1", "FC:DEFB1","p-value:DEFB1",
                       "log2:DEFB4A", "FC:DEFB4A","p-value:DEFB4A",
                       "log2:DEFB103B", "FC:DEFB103B","p-value:DEFB103B",
                       "log2:DEFB104A", "FC:DEFB104A","p-value:DEFB104A",
                       "log2:CAMP", "FC:CAMP","p-value:CAMP",
                       "log2:REG3A", "FC:REG3A","p-value:REG3A",
                       "log2:PLA2G2A", "FC:PLA2G2A","p-value:PLA2G2A",
                       "log2:IL1B", "FC:IL1B","p-value:IL1B",
                       "log2:CXCL8", "FC:CXCL8","p-value:CXCL8",
                       "log2:TNF", "FC:TNF","p-value:TNF" )
#not all genes were found on this platform!
U133PLUS_2_0 <- U133PLUS_2_0[5:nrow(U133PLUS_2_0),]
class(U133PLUS_2_0)
str(U133PLUS_2_0)
U133PLUS_2_0_modified <-data.frame(lapply(U133PLUS_2_0,as.numeric))
U133PLUS_2_0_modified$Perturbations <- U133PLUS_2_0$Perturbations
str(U133A_2_modified)
U133PLUS_2_0_modified

#compounds that match the conditions
AMP_U133PLUS_2_0 <- c("FC.BM2","FC.DEFB1","FC.DEFB4A","FC.DEFB103B","FC.DEFB104A", "FC.CAMP","FC.REG3A","FC.PLA2G2A")
PRO_U133PLUS_2_0 <- c("FC.IL1B","FC.CXCL8","FC.TNF")

results_U133PLUS <- DE_Genevestigator(U133PLUS_2_0_modified,AMP_data = AMP_U133PLUS_2_0,PRO_data = PRO_U133PLUS_2_0)

#what compounds and how many?
selected_words_PLUS <- c()
for (i in 1:length(results_U133PLUS$Perturbations)){
  selected_words_PLUS[i] <- str_extract(results_U133PLUS$Perturbations[i], "([^\\s]+)") 
}
selected_words_PLUS <- selected_words_PLUS[!duplicated(selected_words_PLUS)]
length(selected_words_PLUS)


##U219-4

U219_4 <-read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data/Toxicology-HS_AFFY_U219-4.xlsx")

#modifying the table for analysis 
U219_4 <-U219_4[,-c(1,3,4,5)]
colnames(U219_4) <- c("Perturbations", "log2:BM2", "FC:BM2","p-value:BM2",
                            "log2:DEFB1", "FC:DEFB1","p-value:DEFB1",
                            "log2:CAMP", "FC:CAMP","p-value:CAMP",
                            "log2:REG3A", "FC:REG3A","p-value:REG3A",
                            "log2:PLA2G2A", "FC:PLA2G2A","p-value:PLA2G2A",
                            "log2:IL1B", "FC:IL1B","p-value:IL1B",
                            "log2:CXCL8", "FC:CXCL8","p-value:CXCL8",
                            "log2:TNF", "FC:TNF","p-value:TNF" )
#not all genes were found on this platform!
U219_4 <- U219_4[5:nrow(U219_4),]
class(U219_4)
str(U219_4)
U219_4_modified <-data.frame(lapply(U219_4,as.numeric))
U219_4_modified$Perturbations <- U219_4$Perturbations
str(U219_4_modified)
U219_4_modified

AMP_U219_4 <- c("FC.BM2","FC.DEFB1","FC.CAMP","FC.REG3A","FC.PLA2G2A")
PRO_U219_4 <- c("FC.IL1B","FC.CXCL8","FC.TNF")

View(DE_Genevestigator(U219_4_modified,AMP_data = AMP_U219_4,PRO_data = PRO_U219_4)) # nothing

View(DE_Genevestigator(U219_4_modified,AMP_data = AMP_U219_4,PRO_data = PRO_U219_4,AMP = 0, PRO = 0)) # but works

#nothing!


###LINC-100

LINC_L1000 <-read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data/Toxicology-HS_LINC_L1000-3.xlsx")

#modifying the table for analysis 
LINC_L1000 <-LINC_L1000[,-c(1,3,4,5)]
colnames(LINC_L1000) <- c("Perturbations", 
                          "log2:BM2", "FC:BM2","p-value:BM2",
                      "log2:DEFB1", "FC:DEFB1","p-value:DEFB1",
                      "log2:DEFB4A", "FC:DEFB4A","p-value:DEFB4A",
                      "log2:CAMP", "FC:CAMP","p-value:CAMP",
                      "log2:REG3A", "FC:REG3A","p-value:REG3A",
                      "log2:PLA2G2A", "FC:PLA2G2A","p-value:PLA2G2A",
                      "log2:IL1B", "FC:IL1B","p-value:IL1B",
                      "log2:CXCL8", "FC:CXCL8","p-value:CXCL8",
                      "log2:TNF", "FC:TNF","p-value:TNF" )
#not all genes were found on this platform!
LINC_L1000 <- LINC_L1000[5:nrow(LINC_L1000),]
class(LINC_L1000)
str(LINC_L1000)
LINC_L1000_modified <-data.frame(lapply(LINC_L1000,as.numeric))
LINC_L1000_modified$Perturbations <- LINC_L1000$Perturbations
str(LINC_L1000_modified)
LINC_L1000_modified

AMP_LINC_L1000 <- c("FC.BM2","FC.DEFB1","FC.DEFB4A","FC.CAMP","FC.REG3A","FC.PLA2G2A")
PRO_LINC_L1000 <- c("FC.IL1B","FC.CXCL8","FC.TNF")

results_LINC <-DE_Genevestigator(LINC_L1000_modified,AMP_data = AMP_LINC_L1000,PRO_data = PRO_LINC_L1000)


### maybe let's leave the bit with getting the names for the next step of the analysis, now it'd be good to get the ratios bc at is it the case in LINC_1000 there are way too many compounds!

selected_words <- c()
##ensuring I get the compound only once
for (i in 1:length(results_LINC$Perturbations)){
  selected_words[i] <- str_extract(results_LINC$Perturbations[i], "([^\\s]+)") 
}
selected_words <- selected_words[!duplicated(selected_words)]
length(selected_words) # 270 compounds 


#how to retrieve names and build something out of it?
name_for_storage <- c(paste0("experiment", deparse(substitute(data))))
name_for_storage <-c()

############ Try to represent data for all platforms together, as a heat map. remember one study per compound #######

## 

##### platforms:#####
results_LINC #nrow6809 ncol:28
results_U133A #nrow28 ncol:28
results_U133PLUS #nrow38 ncol:34

# for_index <- results_U133PLUS
# for_index$Index <- c(1:nrow(for_index))
# ##do it together for results_U133A (created below) & results_U133PLUS##
# #created a table in 

Pert_names_before_U133PLUS <-results_U133PLUS$Perturbations #keep pertubrations column, keep names only!

Pert_names_after_U133PLUS <-c()
for (i in 1:length(Pert_names_before_U133PLUS)) {
  Pert_names_after_U133PLUS[i] <- str_extract(Pert_names_before_U133PLUS[i], "([^\\s]+)") 
}

U133PLUS_tobeJoined <-results_U133PLUS[ , grepl( "FC." , names(results_U133PLUS) ) ]

U133PLUS_tobeJoined$Perturbations <- Pert_names_after_U133PLUS

#remove repetitions
U133PLUS_tobeJoined_final <- U133PLUS_tobeJoined[-c(2,3,5,6,7,9,10,15,16,17,22,24,30,31,35),]

colnames(U133PLUS_tobeJoined_final) <- c("BM2", "DEFB1","DEFB4A","DEFB103B","DEFB104A", "CAMP", "REG3A", "PLA2G2A", "IL1B", "IL8", "TNF", "Perturbations")

# done separately for results_U133A

results_U133A$Index <- c(1:nrow(results_U133A))
#remove duplicated studies, when more cells keep cell MCF7, so remove row: 14,17,28,26
results_U133A_for_analysis <- results_U133A[-c(14,17,26,28),]
#remove columns you won't use for the analysis
results_U133A_for_analysis <- results_U133A_for_analysis
Perturbations_to_keep <-results_U133A_for_analysis$Perturbations #keep pertubrations column, keep names only!
for (i in 1:length(Perturbations_to_keep)) {
  Perturbations_to_keep[i] <- str_extract(Perturbations_to_keep[i], "([^\\s]+)") 
}
# perturbation column is lost here!
results_U133A_for_analysis <-results_U133A_for_analysis[ , grepl( "FC." , names(results_U133A_for_analysis ) ) ]

#change names of columns
colnames(results_U133A_for_analysis) <- c("BM2", "DEFB1","DEFB4A", "CAMP", "REG3A", "PLA2G2A", "IL1B", "IL8", "TNF")

results_U133A_for_analysis$Perturbations <-Perturbations_to_keep

#combine the two tables! # they have different dimensions, results_U133A_for_analysis is missing two columns 

U133PLUS_tobeJoined_final
results_U133A_for_analysis
joined_tables <-full_join(U133PLUS_tobeJoined_final, results_U133A_for_analysis) 

#hierarchical analysls- the format of the table is correct, just remove the perturbations

hierarchical_analysis <-joined_tables[,2:(ncol(joined_tables)-1)]
rownames(hierarchical_analysis) <- make.names(joined_tables[,12], unique = TRUE)

unsc_ward <- pheatmap(hierarchical_analysis, clustering_method = "complete", show_colnames = TRUE, show_rownames = TRUE, annotation_legend = FALSE, cluster_cols = FALSE,scale = "row",)

  
# # ## example from data analysis class:
# exprs <- read.csv("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Data analysis/project/data/nanostring.txt", header = TRUE, sep="\t", stringsAsFactors = TRUE)
# dim(exprs)
# str(exprs)
# 
# exprs.num <- as.matrix(exprs[,-(1:2)]) #removing categorical variables
# rownames(exprs.num) <-  make.names(exprs[,1], unique = TRUE)
# exprs.cat <- as.data.frame(exprs[, 2])
# rownames(exprs.cat) <- make.names(exprs[,1], unique = TRUE)

## this hierarchical thing may not be what I need, instead use heatmap
joined_tables

joined_tables_for_visualisation <-gather(data =joined_tables,key = "Gene", value = "Expression_Ratio", -Perturbations)

joined_tables_for_visualisation$Gene <- factor(joined_tables_for_visualisation$Gene, levels = c("BM2", "DEFB1","DEFB4A","DEFB103B", "DEFB104A", "CAMP", "REG3A", "PLA2G2A", "IL1B", "IL8", "TNF"))

ggplot(joined_tables_for_visualisation, aes(x = Gene, y = Perturbations, fill = Expression_Ratio)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis()

##okay now let's work on the LINC data
  #first I need to select only one experiment from each
  #for that I need to extract the concentration, time and cell type, I need to insert them as additional columns


#extract all cell names in LINC:

# bar <- results_LINC$Perturbations[2]
# foo <- "(+)-JQ1 (24h; 0.04uM; A375) / DMSO control (24h; A375)"

LINC_perturbations <-results_LINC$Perturbations
length(LINC_perturbations) #nrow - 6809

LINC_time <- c()
LINC_conc <- c()
LINC_cell <- c()
LINC_compound <- c()

for (i in 1:length(LINC_perturbations)){
  
  linc_perturbations_results <- str_match(as.list(LINC_perturbations[i]), ".+\\((.+)\\;\\s?(.+)\\;\\s?(.+)\\).*\\(.*")
  LINC_time[i] <- linc_perturbations_results[2]
  LINC_conc[i] <- linc_perturbations_results[3]
  LINC_cell[i] <- linc_perturbations_results[4]
  
  LINC_compound[i] <-str_extract(results_LINC$Perturbations[i], "([^\\s]+)") 
  

}

#before you add  new columns, extract only FC columns
results_LINC_1 <-results_LINC[, grepl("FC.", names(results_LINC))]
colnames(results_LINC_1) <- c("BM2","DEFB1", "DEFB4A", "CAMP", "REG3A", "PLA2G2A", "IL1B", "IL8", "TNF")

#now add time, conc, cell type and compound to the result list:
results_LINC_1$Time <-LINC_time
results_LINC_1$Concentration <-LINC_conc
results_LINC_1$Cell <-LINC_cell
results_LINC_1$Perturbations <-LINC_compound

unique(LINC_conc) #"0.04uM" "0.12uM" "0.37uM" "1.11uM" "3.33uM" "10uM"  
unique(LINC_time) # "24h" "3h" 
unique(LINC_cell) #"A375"     "HA1E"     "HT29"     "A549"     "HEPG2"    "MCF7" [7] "HCC515"   "PC3"      "BT20"     "HS578T"   "SKBR3"    "HME1"    "LNCAP"    "MCF10A"   "MDAMB231"

LINC_unique_only<-results_LINC_1[!duplicated(results_LINC_1$Perturbations),]


# hierarchical clustering 
LINC_hierarchy <-LINC_unique_only[,2:9]
rownames(LINC_hierarchy) <- make.names(LINC_unique_only[,13], unique = TRUE)

unsc_ward_LINC <- pheatmap(LINC_hierarchy, clustering_method = "ward.D2", show_colnames = TRUE, show_rownames = FALSE, annotation_legend = FALSE, scale = "row", cluster_cols = FALSE)

# ratios for all tables (unique values for the moment)

U133PLUS_tobeJoined_final
results_U133A_for_analysis
LINC_unique_only

all_platforms_joined <-full_join(joined_tables,LINC_unique_only[,-c(10:12)])

AMP_genes <- c("DEFB1",	"DEFB4A",	"DEFB103B",	"DEFB104A",	"CAMP",	"REG3A",	"PLA2G2A")
PRO_genes <- c("IL1B","IL8","TNF")
all_platforms_joined 
  
ratio_values <- c() 
for (i in 1:nrow(all_platforms_joined)){
  ratio_values[i] <- (sum(all_platforms_joined[i,AMP_genes],na.rm = TRUE)/length(AMP_genes))/(sum(all_platforms_joined[i,PRO_genes],na.rm = TRUE)/length(PRO_genes))
}

all_platforms_joined$Ratios <-ratio_values

ggplot(all_platforms_joined, aes(x = Perturbations, y = Ratios, colour = dplyr::case_when(Ratios > 1 ~ ">1", Ratios < 1 ~ "<1"))) +
  geom_point() + 
  scale_y_continuous(trans = "log10") +
  theme_classic() +
   theme(axis.text.x=element_blank(),
         axis.ticks.x=element_blank())+
  scale_colour_discrete("Ratio") +
  xlab("Compound")+
  ylab("Ratio (AMP/PRO)")
## now I can say, filter the table to only get values with ratio > 1

compounds_with_higher_ratios <-all_platforms_joined %>%
  filter(Ratios > 1) ## 52 compounds!

 #hierarchical clustering of this table?
higher_ratios_hierarchical  <- compounds_with_higher_ratios[,2:(ncol(compounds_with_higher_ratios)-2)]

pheatmap(higher_ratios_hierarchical, clustering_method = "complete", show_colnames = TRUE, show_rownames = FALSE, annotation_legend = FALSE, cluster_cols = FALSE,scale = "row",main = "Compounds with ratio above 1")
      ## so from all of the compounds, these are the ones that have ratio above 1. 

#alternatively what I can do is to get the highest ratios simply and choose the 10-20 first values
sorted_ratios_all_compounds <-all_platforms_joined %>%
  arrange(desc(Ratios))

View(head(sorted_ratios_all_compounds,30))

#################### SEPERATE COMPOUNDS ################
#okay it seems that ibrutinib has the highest ratio so let's look at it in the original table LINC_1000:
dim(LINC_L1000)


LINC_full_table_time <- c()
LINC_full_table_conc <- c()
LINC_full_table_cell <- c()
LINC_full_table_compound <- c()

for (i in 1:length(LINC_L1000$Perturbations)){
  
  linc_perturbations_results <- str_match(as.list(LINC_L1000$Perturbations[i]), ".+\\((.+)\\;\\s?(.+)\\;\\s?(.+)\\).*\\(.*")
  LINC_full_table_time[i] <- linc_perturbations_results[2]
  LINC_full_table_conc[i] <- linc_perturbations_results[3]
  LINC_full_table_cell[i] <- linc_perturbations_results[4]
  
  LINC_full_table_compound[i] <-str_extract(LINC_L1000$Perturbations[i], "([^\\s]+)") 
  
  
}

LINC_L1000_FC <- LINC_L1000[, grepl("FC.",names(LINC_L1000))]
LINC_L1000_FC <-data.frame(lapply(LINC_L1000_FC,as.numeric))
str(LINC_L1000_FC)

LINC_L1000_FC$Compound <- LINC_full_table_compound
LINC_L1000_FC$Time <-LINC_full_table_time
LINC_L1000_FC$Concentration <- LINC_full_table_conc
LINC_L1000_FC$Cell <-LINC_full_table_cell

ibrutinib <- LINC_L1000_FC %>%
  filter(Compound == "ibrutinib")

unique(ibrutinib$Concentration) #"0.04uM" "0.12uM" "0.37uM" "1.11uM" "3.33uM" "10uM" 
unique(ibrutinib$Cell) #"A375"   "A549"   "HA1E"   "HCC515" "HEPG2"  "HT29"   "MCF7"   "PC3" 
# all cell types tested for each concentration, let's start off with A375 maybe?

ibrutinib_gathered <- gather(data = ibrutinib,key = "Gene", value = "Expression_Ratio", -c(Compound, Time, Concentration,Cell))

ibrutinib_gathered$Gene <- factor(ibrutinib_gathered$Gene, levels = c("FC.BM2", "FC.DEFB1","FC.DEFB4A","FC.CAMP", "FC.REG3A", "FC.PLA2G2A", "FC.IL1B", "FC.CXCL8", "FC.TNF"))

ibrutinib_gathered$Concentration <- factor(ibrutinib_gathered$Concentration, levels = c("0.04uM", "0.12uM", "0.37uM", "1.11uM", "3.33uM", "10uM"))

A375 <-ggplot(ibrutinib_gathered[ibrutinib_gathered$Cell == "A375",],mapping = aes(x = Gene, y = Expression_Ratio, fill =Concentration)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() 

A549 <-ggplot(ibrutinib_gathered[ibrutinib_gathered$Cell == "A549",],mapping = aes(x = Gene, y = Expression_Ratio, fill =Concentration)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() 

ggplot(ibrutinib_gathered[ibrutinib_gathered$Cell == "HA1E",],mapping = aes(x = Gene, y = Expression_Ratio, fill =Concentration)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() 

ggplot(ibrutinib_gathered[ibrutinib_gathered$Cell == "HCC515",],mapping = aes(x = Gene, y = Expression_Ratio, fill =Concentration)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() 

ggplot(ibrutinib_gathered[ibrutinib_gathered$Cell == "HEPG2",],mapping = aes(x = Gene, y = Expression_Ratio, fill =Concentration)) +
  geom_bar(stat = "identity",position = "dodge") +
  theme_classic() +
  coord_flip()

#line plot, confusing!
ggplot(ibrutinib_gathered[ibrutinib_gathered$Gene == "FC.DEFB4A",],mapping = aes(x = Cell, y = Expression_Ratio, colour =Concentration)) +
 geom_point() +
  geom_line(aes(group = Concentration)) +
theme_classic() +
  ylab("HBD2:Expression Ratio")

### this is a niec one! ti shows you there are differences among cell types
ggplot(ibrutinib_gathered[ibrutinib_gathered$Gene == "FC.CAMP",],mapping = aes(x = Cell, y = Expression_Ratio)) +
  geom_boxplot() +
  geom_jitter(aes(colour = Concentration)) +
  theme_classic() +
  ylab("REG3A:Expression Ratio")

#Condtions: only 24h,different types of concentrations, different cells

## can also represent in the PCA way!




