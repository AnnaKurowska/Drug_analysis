---
title: "Genestigator-based analysis"
author: "Anna Kurowska"
date: "21/05/2022"
output: ioslides_presentation
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
```


#1.Loading libraries
```{r, echo= = FALSE}
library(rlang)
library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(pheatmap)
library(ggpubr)
library(hash) # good question
```

#2.Loading and processing datasets from 4 platforms

In the following section, download the original datasets from U133A, U133APLUS, U219, LINCS and prepare them for the analaysis.

Platform U133A:
```{r}
#load the dataset
U133A <- read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/Downloaded_tables/HS_AFFY_U133A.xlsx")

#remove unnecessary columns
U133A <-U133A[,-c(1,3,4,5)]
#change the column names for easier handling
colnames(U133A) <- c("Perturbations", "FC:BM2","p-value:BM2",
                       "FC:HBD1","p-value:HBD1",  
                       "FC:HBD2","p-value:HBD2",
                       "FC:LL37","p-value:LL37",
                       "FC:REG3A","p-value:REG3A",
                       "FC:PLA2G2A","p-value:PLA2G2A",
                       "FC:PGLYRP1","p-value:PGLYRP1",
                       "FC:PGLYRP4","p-value:PGLYRP4",
                       "FC:IL1B","p-value:IL1B",
                       "FC:IL8","p-value:IL8",
                       "FC:CCL20","p-value:CCL20",
                       "FC:TNF","p-value:TNF" )

#remove unnecessary rows
U133A <- U133A[4:nrow(U133A),]
#numbers values are in the character format, convert character values into numeric values
U133A_modified <-data.frame(lapply(U133A,as.numeric))
#since conversion into numbers will remove the character Perturbations columns that holds the name of the experiments, add it manually to a new table
U133A_modified$Perturbations <- U133A$Perturbations
#check  dimensions
dim(U133A_modified)

#from the experiments column use REGEX to only retrieve the names of the molecule used in the experiment. Both full names of the experiments and molecule names will be kept
Pert_names_before_U133 <-U133A_modified$Perturbations

Pert_names_after_U133<-c()
for (i in 1:length(Pert_names_before_U133)) {
  Pert_names_after_U133[i] <- str_extract(Pert_names_before_U133[i], "([^\\s]+)") 
}
#select columns that contain FC (fold change) values
U133A_modified <-U133A_modified[ , grepl( "FC." , names(U133A_modified) ) ]
#bind experiments and molecule columns
U133A_modified <-cbind(as.data.frame(Pert_names_after_U133),as.data.frame(Pert_names_before_U133), U133A_modified)
#change the final naming of the columns
colnames(U133A_modified) <- c("Compounds", "Original","BM2", "HBD1","HBD2", "LL37", "REG3A", "PLA2G2A","PGLYRP1","PGLYRP4", "IL1B", "IL8", "CCL20", "TNF")

```

Platform U133APLUS
```{r}
#load the dataset
U133A_PLUS <- read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/Downloaded_tables/HS_AFFY_U133PLUS_2.xlsx")
#remove unnecessary columns
U133A_PLUS <-U133A_PLUS[,-c(1,3,4,5)]
#change the column names for easier handling
colnames(U133A_PLUS) <- c("Perturbations", 
                       "FC:BM2","p-value:BM2",
                       "FC:HBD1","p-value:HBD1",  
                       "FC:HBD2","p-value:HBD2",  
                       "FC:HBD3","p-value:HBD3",
                       "FC:HBD4","p-value:HBD4",
                       "FC:LL37","p-value:LL37",
                       "FC:REG3A","p-value:REG3A",
                       "FC:PLA2G2A","p-value:PLA2G2A",
                       "FC:PGLYRP1","p-value:PGLYRP1",
                       "FC:PGLYRP2","p-value:PGLYRP2",
                       "FC:PGLYRP3","p-value:PGLYRP3",
                       "FC:PGLYRP4","p-value:PGLYRP4",
                       "FC:IL1B","p-value:IL1B",
                       "FC:IL8","p-value:IL8",
                       "FC:CCL20","p-value:CCL20",
                       "FC:TNF","p-value:TNF" )
#remove unnecessary rows
U133A_PLUS <- U133A_PLUS[4:nrow(U133A_PLUS),]
#numbers values are in the character format, convert character values into numeric values
U133APLUS_modified <-data.frame(lapply(U133A_PLUS,as.numeric))
#since conversion into numbers will remove the character Perturbations columns that holds the name of the experiments, add it manually to a new table
U133APLUS_modified$Perturbations <- U133A_PLUS$Perturbations
#check  dimensions
dim(U133APLUS_modified)

#from the experiments column use REGEX to only retrieve the names of the molecule used in the experiment. Both full names of the experiments and molecule names will be kept
Pert_names_before_U133PLUS <-U133APLUS_modified$Perturbations 

Pert_names_after_U133PLUS<-c()
for (i in 1:length(Pert_names_before_U133PLUS)) {
  Pert_names_after_U133PLUS[i] <- str_extract(Pert_names_before_U133PLUS[i], "([^\\s]+)") 
}
#select columns that contain FC values
U133APLUS_modified <-U133APLUS_modified[, grepl("FC.",names(U133APLUS_modified)) ]

#bind experiments and molecule columns
U133APLUS_modified <-cbind(as.data.frame(Pert_names_after_U133PLUS),as.data.frame(Pert_names_before_U133PLUS), U133APLUS_modified)

#change the final naming of the columns
colnames(U133APLUS_modified) <- c("Compounds", "Original","BM2", "HBD1","HBD2","HBD3","HBD4", "LL37", "REG3A", "PLA2G2A","PGLYRP1","PGLYRP2","PGLYRP3","PGLYRP4", "IL1B", "IL8", "CCL20", "TNF")

```

Platform U219
```{r}
#load the dataset
U219 <- read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/Downloaded_tables/HS_AFFY_U219.xlsx")
#remove unnecessary columns
U219 <-U219[,-c(1,3,4,5)]
#change the column names for easier handling
colnames(U219) <- c("Perturbations", 
                       "FC:BM2","p-value:BM2",
                       "FC:HBD1","p-value:HBD1",  
                       "FC:LL37","p-value:LL37",
                       "FC:REG3A","p-value:REG3A",
                       "FC:PLA2G2A","p-value:PLA2G2A",
                       "FC:PGLYRP1","p-value:PGLYRP1",
                       "FC:PGLYRP2","p-value:PGLYRP2",
                       "FC:PGLYRP3","p-value:PGLYRP3",
                       "FC:PGLYRP4","p-value:PGLYRP4",
                       "FC:IL1B","p-value:IL1B",
                       "FC:IL8","p-value:IL8",
                       "FC:CCL20","p-value:CCL20",
                       "FC:TNF","p-value:TNF" )

#remove unnecessary rows
U219 <- U219[4:nrow(U219),]
#numbers values are in the character format, convert character values into numeric values
U219_modified <-data.frame(lapply(U219,as.numeric))

#since conversion into numbers will remove the character Perturbations columns that holds the name of the experiments, add it manually to a new table
U219_modified$Perturbations <- U219$Perturbations
#check  dimensions
dim(U219_modified)

#from the experiments column use REGEX to only retrieve the names of the molecule used in the experiment. Both full names of the experiments and molecule names will be kept
Pert_names_before_U219 <-U219_modified$Perturbations 

Pert_names_after_U219<-c()
for (i in 1:length(Pert_names_before_U219)) {
  Pert_names_after_U219[i] <- str_extract(Pert_names_before_U219[i], "([^\\s]+)") 
}

#select columns that contain FC values
U219_modified <-U219_modified[, grepl("FC.",names(U219_modified)) ]

#bind experiments and molecule columns
U219_modified <-cbind(as.data.frame(Pert_names_after_U219),as.data.frame(Pert_names_before_U219), U219_modified)

#change the final naming of the columns
colnames(U219_modified) <- c("Compounds", "Original","BM2", "HBD1", "LL37", "REG3A", "PLA2G2A","PGLYRP1","PGLYRP2","PGLYRP3","PGLYRP4", "IL1B", "IL8", "CCL20", "TNF")
```

Platform LINCS
```{r}
#load the dataset
LINC <- read_excel("/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/Downloaded_tables/HS_LINC_L1000.xlsx")

#remove unnecessary columns
LINC <-LINC[,-c(1,3,4,5)]
#change the column names for easier handling
colnames(LINC) <- c("Perturbations", 
                       "FC:BM2","p-value:BM2",
                       "FC:HBD1","p-value:HBD1", 
                       "FC:HBD2","p-value:HBD2", 
                       "FC:LL37","p-value:LL37",
                       "FC:REG3A","p-value:REG3A",
                       "FC:PLA2G2A","p-value:PLA2G2A",
                       "FC:PGLYRP1","p-value:PGLYRP1",
                       "FC:PGLYRP4","p-value:PGLYRP4",
                       "FC:IL1B","p-value:IL1B",
                       "FC:IL8","p-value:IL8",
                       "FC:CCL20","p-value:CCL20",
                       "FC:TNF","p-value:TNF" )
#remove unnecessary rows
LINC <- LINC[4:nrow(LINC),]
#numbers values are in the character format, convert character values into numeric values
LINC_modified <-data.frame(lapply(LINC,as.numeric))

#since conversion into numbers will remove the character Perturbations columns that holds the name of the experiments, add it manually to a new table
LINC_modified$Perturbations <- LINC$Perturbations

#check  dimensions
dim(LINC_modified)

#from the experiments column use REGEX to only retrieve the names of the molecule used in the experiment. Both full names of the experiments and molecule names will be kept
Pert_names_before_LINC <-LINC_modified$Perturbations 

Pert_names_after_LINC<-c()
for (i in 1:length(Pert_names_before_LINC)) {
  Pert_names_after_LINC[i] <- str_extract(Pert_names_before_LINC[i], "([^\\s]+)") 
}

#select columns that contain FC values
LINC_modified <-LINC_modified[, grepl("FC.",names(LINC_modified)) ]

#bind experiments and molecule columns
LINC_modified <-cbind(as.data.frame(Pert_names_after_LINC),as.data.frame(Pert_names_before_LINC), LINC_modified)

#change the final naming of the columns
colnames(LINC_modified) <- c("Compounds", "Original","BM2", "HBD1","HBD2", "LL37", "REG3A", "PLA2G2A","PGLYRP1","PGLYRP4", "IL1B", "IL8", "CCL20", "TNF")

```

#3. Find the compounds inducing the highest gene expression for each dataset 

Separate dataset for each platform into separate tables for each gene, do the analysis of each gene to identify compounds.

Function that selects molecules that upregulate the expression of at least one AMP in at least one experiment. Positive molecules are selected for further analysis. Save the output into excel tables that are then inspected and annotated manually.
```{r}
find_expression <- function(data, gene){
  
  ordered_data <- data %>%
  dplyr::filter(!!rlang::parse_expr(gene) >= 2) %>% 
  dplyr::select(c(Compounds, Original, gene, IL1B,IL8, CCL20,TNF)) %>%
  arrange(desc(!!rlang::parse_expr(gene)))

  return(ordered_data)
} 
```

Platform U133A
```{r}
#HBD1

U133A_HBD1 <-find_expression(U133A_modified,"HBD1")
  #none

#HBD2

U133A_HBD2 <- find_expression(U133A_modified,"HBD2")

#write.csv(U133A_HBD2,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133A_HBD2.csv")

#LL37

U133A_LL37 <-find_expression(U133A_modified,"LL37")

#write.csv(U133A_LL37,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133A_LL37.csv")

#REG3A

U133A_REG3A <- find_expression(U133A_modified,"REG3A")
  #none

#PLA2G2A

U133A_PLA2G2A <- find_expression(U133A_modified,"PLA2G2A")

#write.csv(U133A_PLA2G2A,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133A_PLA2G2A.csv")

#PGLYRP1

U133A_PGLYRP1 <- find_expression(U133A_modified,"PGLYRP1")

#write.csv(U133A_PGLYRP1,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133A_PGLYRP1.csv")

#PGLYRP4

U133A_PGLYRP4 <- find_expression(U133A_modified,"PGLYRP4")

#write.csv(U133A_PGLYRP4,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133A_PGLYRP4.csv")
```

Platform U133APLUS
```{r}

#"HBD1

U133APLUS_HBD1 <- find_expression(U133APLUS_modified,"HBD1")

#write.csv(U133APLUS_HBD1,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_HBD1.csv")

#HBD2

U133APLUS_HBD2 <- find_expression(U133APLUS_modified,"HBD2")

#write.csv(U133APLUS_HBD2,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_HBD2.csv")

#HBD3

U133APLUS_HBD3 <- find_expression(U133APLUS_modified,"HBD3")

#write.csv(U133APLUS_HBD3,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_HBD3.csv")

#HBD4

U133APLUS_HBD4 <- find_expression(U133APLUS_modified,"HBD4")

#write.csv(U133APLUS_HBD4,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_HBD4.csv")

#LL37

U133APLUS_LL37 <- find_expression(U133APLUS_modified,"LL37")

#write.csv(U133APLUS_LL37,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_LL37.csv")

#REG3A

U133APLUS_REG3A <- find_expression(U133APLUS_modified,"REG3A")

#write.csv(U133APLUS_REG3A,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_REG3A.csv")

#PLA2G2A 

U133APLUS_PLA2G2A <- find_expression(U133APLUS_modified,"PLA2G2A")

#write.csv(U133APLUS_PLA2G2A,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_PLA2G2A.csv")

#PGLYRP1

U133APLUS_PGLYRP1 <-find_expression(U133APLUS_modified,"PGLYRP1")
  #none

#PGLYRP2

U133APLUS_PGLYRP2 <-find_expression(U133APLUS_modified,"PGLYRP2")

#write.csv(U133APLUS_PGLYRP2,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_PGLYRP2.csv")

#PGLYRP3

U133APLUS_PGLYRP3 <-find_expression(U133APLUS_modified,"PGLYRP3")
  #none

#PGLYRP4

U133APLUS_PGLYRP4 <-  find_expression(U133APLUS_modified,"PGLYRP4")

#write.csv(U133APLUS_PGLYRP4,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/U133APLUS_PGLYRP4.csv")

```

Platform U219
```{r}

#HBD1
find_expression(U219_modified, "HBD1")
  #none

#LL37
find_expression(U219_modified, "LL37")
  #none

#REG3A
find_expression(U219_modified, "REG3A")
  #none

#PLA2G2A
find_expression(U219_modified, "PLA2G2A")
  #none

#PGLYRP1
find_expression(U219_modified, "PGLYRP1")
  #none

#PGLYRP2
find_expression(U219_modified, "PGLYRP2")
  #none

#PGLYRP3
find_expression(U219_modified, "PGLYRP3")
  #none

#PGLYRP4
find_expression(U219_modified, "PGLYRP4")
  #none
```

LINC_modified
```{r}

#HBD1

LINC_HBD1 <- find_expression(LINC_modified, "HBD1")

#write.csv(LINC_HBD1,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/LINC_HBD1.csv")

unique_LINC_HBD1 <- LINC_HBD1[!duplicated(LINC_HBD1$Compounds),]

#write.csv(unique_LINC_HBD1,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/unique_LINC_HBD1.csv")

#HBD2

LINC_HBD2 <- find_expression(LINC_modified, "HBD2")

#write.csv(LINC_HBD2,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/LINC_HBD2.csv")

unique_LINC_HBD2 <- LINC_HBD2[!duplicated(LINC_HBD2$Compounds),]

#write.csv(unique_LINC_HBD2,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/unique_LINC_HBD2.csv")

#LL37

LINC_LL37 <- find_expression(LINC_modified, "LL37")

#write.csv(LINC_LL37,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/LINC_LL37.csv")

unique_LINC_LL37 <- LINC_LL37[!duplicated(LINC_LL37$Compounds),]

#write.csv(unique_LINC_LL37,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/unique_LINC_LL37.csv")

#REG3A

LINC_REG3A <- find_expression(LINC_modified, "REG3A")

#write.csv(LINC_REG3A,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/LINC_REG3A.csv")

unique_LINC_REG3A <- LINC_REG3A[!duplicated(LINC_REG3A$Compounds),]

#write.csv(unique_LINC_REG3A,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/unique_LINC_REG3A.csv")

#PLA2G2A

LINC_PLA2G2A <- find_expression(LINC_modified, "PLA2G2A")

#write.csv(LINC_PLA2G2A,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/LINC_PLA2G2A.csv")

unique_LINC_PLA2G2A <- LINC_PLA2G2A[!duplicated(LINC_PLA2G2A$Compounds),]

#write.csv(unique_LINC_PLA2G2A,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/unique_LINC_PLA2G2A.csv")

#PGLYRP1

LINC_PGLYRP1 <- find_expression(LINC_modified, "PGLYRP1")

#write.csv(LINC_PGLYRP1,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/LINC_PGLYRP1.csv")

unique_LINC_PGLYRP1 <- LINC_PGLYRP1[!duplicated(LINC_PGLYRP1$Compounds),]

#write.csv(unique_LINC_PGLYRP1,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/unique_LINC_PGLYRP1.csv")

#PGLYRP4

LINC_PGLYRP4 <- find_expression(LINC_modified, "PGLYRP4")

#write.csv(LINC_PGLYRP4,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/LINC_PGLYRP4.csv")

unique_LINC_PGLYRP4 <- LINC_PGLYRP4[!duplicated(LINC_PGLYRP4$Compounds),]

#write.csv(unique_LINC_PGLYRP4,"/Users/Ania/Desktop/Szkoła/ENS /2nd_year/Stage/Practical/Genevestigator/Data_17_02_22_current/output_tables_for_excel/unique_LINC_PGLYRP4.csv")


```

#4. Cell types in each platform: general description of the data

```{r}

#U133A
U13AA_cells <- c()

for (i in 1:nrow(U133A_modified)){
  U13AA_cells[i] <- str_match(U133A_modified[i,2],".*treated (.*) cell sample")[[2]]
}
#mode: str_match("cobalt chloride study 1 (100uM) / vehicle (DMSO) treated MCF7 cell sample",".*treated (.*) cell sample")[[2]]

U13AA_cells <-as.data.frame(table(U13AA_cells))
colnames(U13AA_cells) <- c("Cell", "Freq")

#U133APLUS 

U13AAPLUS_cells<- c() ##This retrieves more, maybe just use that?#################

for (i in 1:nrow(U133APLUS_modified)){
  U13AAPLUS_cells[i] <- str_match(U133APLUS_modified[i,2],".*treated (.*) sample")[[2]]
}

 #remove NA
  U13AAPLUS_cells <- U13AAPLUS_cells[!is.na(U13AAPLUS_cells)]
   #numer of rows retrieved here
  length(U13AAPLUS_cells) 

U13AAPLUS_cells <-as.data.frame(table(U13AAPLUS_cells))
colnames(U13AAPLUS_cells) <- c("Cell", "Freq")
  
#U219

U219_cells <- c()

for (i in 1:nrow(U219_modified)){
  U219_cells[i] <- str_match(U219_modified[i,2],".*\\(.*; (.*)\\) \\/ .*")[[2]]
}

#model: str_match("trenbolone study 1 (100nM; Ishikawa) / vehicle (DMSO) treated Ishikawa cell sample",".*\\(.*; (.*)\\) \\/ .*")[[2]]

U219_cells <-as.data.frame(table(U219_cells))
colnames(U219_cells) <- c("Cell", "Freq")


# LINCS
LINC_cell <- as.data.frame(table(LINC_cell))
colnames(LINC_cell) <- c("Cell", "Freq")

cell_types <-rbind(U13AA_cells,U13AAPLUS_cells,U219_cells, LINC_cell)

cell_types %>%
  group_by(Cell) %>%
  summarise(Frequency = sum(Freq)) 

##in total, 24441 experiments

#HT-29
1611/24441*100  6.591383

```

#5.Analysis of LINC compounds in HT-29 cell type. 

Creating a summary table for HT-29 cells treated with HT-29 cells. The goal is to calculate the frequency of each gene upregulation (FC >= 2) by the molecules from the Genevestigator database. Since only the LINCS platform contained experiments done in colonic cell lines (HT-29) only data from this platform is analysed.

```{r}
# Set up empty variables
LINC_perturbations <- LINC_modified$Original
  
LINC_time <- c()
LINC_conc <- c()
LINC_cell <- c()
LINC_compound <- c()

#for every molecule acquired from the LINCS platform retrieve the informations regarding the time, cell, and concentration used. Possible to get all of this information as the experiments are annotated the same way.
for (i in 1:length(LINC_perturbations)){

  linc_perturbations_results <- str_match(as.list(LINC_perturbations[i]), ".+\\((.+)\\;\\s?(.+)\\;\\s?(.+)\\).*\\(.*")
  LINC_time[i] <- linc_perturbations_results[2]
  LINC_conc[i] <- linc_perturbations_results[3]
  LINC_cell[i] <- linc_perturbations_results[4]
  LINC_compound[i] <-str_extract(LINC_perturbations[i], "([^\\s]+)") 
  
}
#bind all the variables together with the LINCS dataframe 
LINC_compounds <-cbind(as.data.frame(LINC_cell),as.data.frame(LINC_conc),as.data.frame(LINC_time), as.data.frame(LINC_modified))

#create HT-29 dataframe only
LINC_compounds_HT29 <-LINC_compounds %>%
  filter(LINC_cell == "HT29") 

##Quantify the frequency of each AMP gene upregulation by each molecule
#Put the data in the narrow format
LINC_HT29_f_gathered <- gather(LINC_compounds_HT29,key = Gene, value = Expression_ratio, -c(1:5)) 

AMP_LINC <- c("HBD1", "HBD2", "LL37", "REG3A", "PLA2G2A", "PGLYRP1", "PGLYRP4")

HT29_freq_each_gene <-LINC_HT29_f_gathered %>%
  filter(Gene %in% AMP_LINC) %>%
  group_by(Compounds,Gene) %>%
  summarise(Gene_Frequency = sum((Expression_ratio >= 2) == TRUE)) %>%
  pivot_wider(names_from = Gene,values_from = Gene_Frequency)
```

Extra: written summary of the AMP genes upregulated by molecules tested on HT-29 cells
```{r}
#for each molecule identify whether at least one AMP gene is upergulated
HT29_logic <- as.data.frame(HT29_freq_each_gene[,AMP_LINC] >= 1)

#for each molecule retrieve names of genes it upregultes
HT29_genes_concatenated <- list()
for (i in 1:nrow(HT29_logic)) {
  numberss <- which(HT29_logic[i,] == TRUE)
  HT29_genes_concatenated[[i]] <- AMP_LINC[numberss]
}

#we have a list of list, name each list containing names of AMP genes with the molecule it corresponds to and convert it into a dataframe
names(HT29_genes_concatenated) <- HT29_freq_each_gene$Compounds
nm_HT29 <- names(HT29_genes_concatenated)
#unlist the list
result_HT29 <- lapply(unique(nm_HT29), function(n) unname(unlist(HT29_genes_concatenated[nm_HT29 %in% n])))
#if a gene was upregulated multiple times
names(result_HT29) <- unique(nm_HT29)
result_HT29 <-lapply(result_HT29, function(x) toString(unique(x)))
result_HT29 <-as.data.frame(unlist(result_HT29))

```

#6.Check the frequency of experiment upregulation

It is of interest to understand how many molecules upregualated at least one gene in the experiments. 

Platform U133A
```{r} 
#U133A

up_experiments_U133A <- c()

for (i in 1:nrow(U133A_modified)){
  for (j in 3:ncol(U133A_modified)) {
  if (U133A_modified[i,j] >= 2 || U133A_modified[i,j] <= -2){
    up_experiments_U133A[i] <- U133A_modified[i,2]
    }
  }
}

up_experiments_U133A <- up_experiments_U133A[!is.na(up_experiments_U133A)]
#calculating proportion of experiments that upregulate at least one gene
length(up_experiments_U133A)/nrow(U133A_modified)* 100
# 245 out of 1645 (14.88457%)

```

Platform U133APLUS
```{r}
#U133APLUS
up_experiments_U133APLUS <- c()

for (i in 1:nrow(U133APLUS_modified)){
  for (j in 3:ncol(U133APLUS_modified)) {
  if (U133APLUS_modified[i,j] >= 2 || U133APLUS_modified[i,j] <= -2){
    up_experiments_U133APLUS[i] <- U133APLUS_modified[i,2]
    }
  }
}

up_experiments_U133APLUS <- up_experiments_U133APLUS[!is.na(up_experiments_U133APLUS)]
#calculating proportion of experiments that upregulate at least one gene
length(up_experiments_U133APLUS)/nrow(U133APLUS_modified)* 100
# 239 out of 1107 (22.49322%)
```

Platform U219
```{r}
#U219
up_experiments_U219 <- c()

for (i in 1:nrow(U219_modified)){
  for (j in 3:ncol(U219_modified)) {
  if (U219_modified[i,j] >= 2 || U219_modified[i,j] <= -2){
    up_experiments_U219[i] <- U219_modified[i,2]
    }
  }
}

up_experiments_U219 <- up_experiments_U219[!is.na(up_experiments_U219)]
#calculating proportion of experiments that upregulate at least one gene
length(up_experiments_U219)/nrow(U219_modified)* 100

```

Platform LINCS
```{r}
#LINCS
up_experiments_LINCS <- c()

for (i in 1:nrow(LINC_modified)){
  for (j in 3:ncol(LINC_modified)) {
  if (LINC_modified[i,j] >= 2 || LINC_modified[i,j] <= -2){
    up_experiments_LINCS[i] <- LINC_modified[i,2]
    }
  }
}

up_experiments_LINCS <- up_experiments_LINCS[!is.na(up_experiments_LINCS)]

#calculating proportion of experiments that upregulate at least one gene
length(up_experiments_LINCS)/nrow(LINC_modified)* 100
# 20288 out of 21413 (94.74618%)
```

#7. Calculating max and min DE values induced by molecules in each platform 


Platform LINCS
```{r}
LINC_max <- as.data.frame(apply(LINC_modified[,3:ncol(LINC_modified)],MARGIN = 2, max))
LINC_max["HBD3",] <- NA #add rows of genes that are not measured in this platform, otherwise you can't join data frames
LINC_max["HBD4",] <- NA
LINC_max["PGLYRP2",] <- NA
LINC_max["PGLYRP3",] <- NA
LINC_min <- as.data.frame(apply(LINC_modified[,3:ncol(LINC_modified)],MARGIN = 2, min))
LINC_min["HBD3",] <- NA
LINC_min["HBD4",] <- NA
LINC_min["PGLYRP2",] <- NA
LINC_min["PGLYRP3",] <- NA
```

Platform U133A
```{r}
U133A_max <- as.data.frame(apply(U133A_modified[,3:ncol(U133A_modified)],MARGIN = 2, max))
U133A_max["HBD3",] <- NA #add rows of genes that are not measured in this platform, otherwise you can't join data frames
U133A_max["HBD4",] <- NA
U133A_max["PGLYRP2",] <- NA
U133A_max["PGLYRP3",] <- NA
U133A_min <- as.data.frame(apply(U133A_modified[,3:ncol(U133A_modified)],MARGIN = 2, min))
U133A_min["HBD3",] <- NA
U133A_min["HBD4",] <- NA
U133A_min["PGLYRP2",] <- NA
U133A_min["PGLYRP3",] <- NA
```


Platform U133APLUS
```{r}
U133APLUS_max <- as.data.frame(apply(U133APLUS_modified[,3:ncol(U133APLUS_modified)],MARGIN = 2, max))
U133APLUS_min <- as.data.frame(apply(U133APLUS_modified[,3:ncol(U133APLUS_modified)],MARGIN = 2, min))
```

Platform U219
```{r}
U129_max <- as.data.frame(apply(U219_modified[,3:ncol(U219_modified)],MARGIN = 2, max))
U129_max["HBD2",] <- NA #add rows of genes that are not measured in this platform, otherwise you can't join data frames
U129_max["HBD3",] <- NA
U129_max["HBD4",] <- NA
U129_min <- as.data.frame(apply(U219_modified[,3:ncol(U219_modified)],MARGIN = 2, min))
U129_min["HBD2",] <- NA
U129_min["HBD3",] <- NA
U129_min["HBD4",] <- NA

```

Join platforms
```{r}
platforms_mix_man <- cbind(U133A_max, 
                          U133A_min,
                          U133APLUS_max,
                          U133APLUS_min,
                          U129_max,
                          U129_min,
                          LINC_max,
                          LINC_min)

colnames(platforms_mix_man) <- c("U133A_max","U133A_min", "U133APLUS_max", "U133APLUS_min", "U219_max","U219_min","LINC_max", "LINC_min" )

platforms_mix_man <- platforms_mix_man[c("BM2", "HBD1", "HBD2", "HBD3", "HBD4", "LL37","REG3A", "PLA2G2A","PGLYRP1", "PGLYRP2", "PGLYRP3", "PGLYRP4", "IL8", "IL1B", "CCL20", "TNF"),]

```

Prepare the table for visualisation
```{r}
platforms_mix_man <- rownames_to_column(platforms_mix_man, var = "Gene") #put rownames into a column

platforms_mix_man_gathered <- gather(platforms_mix_man, key = Category, value = Value,-Gene) #create long data

#create new 2 columns: 1 with platform, 2 with max/min group so that they can be differentiated in the graph
platforms <- c()  
for (i in 1:nrow(platforms_mix_man_gathered)){
  platforms[i] <- str_split(platforms_mix_man_gathered[i,"Category"],pattern = "_")[[1]][1]

}
max_min_group <- c()  
for (i in 1:nrow(platforms_mix_man_gathered)){
  max_min_group[i] <- str_split(platforms_mix_man_gathered[i,"Category"],pattern = "_")[[1]][2]

}
platforms_mix_man_gathered$Platform <- platforms
platforms_mix_man_gathered$Group <- max_min_group

#plot in ggplot2
ggplot(platforms_mix_man_gathered, aes(x = Gene, y = Value, colour = Platform)) +
  geom_point() +
  coord_flip() +
  theme_classic() +
  geom_line(aes(group = Platform, colour = Platform))

value1 <- abs(rnorm(26))*2
data <- data.frame(
  x=LETTERS[1:26], 
  value1=value1, 
  value2=value1+1+rnorm(26, sd=1) 
)

data <- data %>% 
  rowwise() %>% 
  mutate( mymean = mean(c(value1,value2) )) %>% 
  arrange(mymean) %>% 
  mutate(x=factor(x, x))

```
