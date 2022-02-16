library(webchem)
library(caTools)
library(BMS)
library(dplyr)
library(tidyverse)
library(devtools)
install_github("CDK-R/rpubchem", dependencies=TRUE)
library(rpubchem)
library(fingerprint)
library(tsne)
library(ggplot2)

lc50
lc50_sub <- lc50[1:15, ]
lc50_sub$inchikey <- cts_convert(lc50_sub$cas, from = "CAS", to = "InChIKey", choices = 1, verbose = FALSE)
lc50_sub$inchikey <- as.character(lc50_sub$inchikey)
x <- get_cid(lc50_sub$inchikey, from = "inchikey", match = "first", verbose = FALSE)
lc50_sub2 <- full_join(lc50_sub, x, by = c("inchikey" = "query"))

y <- pc_prop(lc50_sub2$cid, properties = c("IUPACName", "XLogP", "Fingerprint2D"))
y$CID <- as.character(y$CID)
lc50_sub3 <- full_join(lc50_sub2, y, by = c("cid" = "CID"))

# finger <-lc50_sub3$Fingerprint2D[1]
# finger_1 <- base64decode(finger, what = "character")
# finger_2 <- hex2bin(enc2utf8(finger_1))
# result <- hex2bin(as.character(base64_dec(finger)))
# binary <- c()
# for (i in 1:nrow(lc50_sub3)) {
#   binary[i] <- base64decode(lc50_sub3$Fingerprint2D[i], what = "character")
#   binary[i] <-hex2bin(as.character(binary[i]))
# }
# what <- decodeCACTVS(finger)  ## now it can be manipulated using the fingerprint package
# what_is_thats <- fp.factor.matrix(list(what))  ## converts fingerprint into factors 1 and 0
# as.character(what_is_thats)

binary <- list()
for (i in 1:nrow(lc50_sub3)) {
  binary[i] <- as.character(decodeCACTVS(lc50_sub3$Fingerprint2D[i]))
}
###turning the character string into a list of integers
matrix_new <-matrix(0,nrow = 15, ncol = 881)
for (i in 1:length(binary)){
  matrix_new[i,] <- as.numeric(t(as.data.frame(strsplit(binary[[i]], ""))))
}

CIDs <- lc50_sub3$cid
compound_data <-as.data.frame(matrix_new)
#compound_data$IDs <-CIDs  not useful now
compound_data$Compound = c("Compound1", "Compound2","Compound3","Compound4","Compound5","Compound6","Compound7", "Compound8","Compound9","Compound10", "Compound11","Compound12", "Compound13", "Compound14", "Compound15")

## random expression levels 
diff_expression <- data.frame(
  Compound = c("Compound1", "Compound2","Compound3","Compound4","Compound5","Compound6","Compound7", "Compound8","Compound9","Compound10", "Compound11","Compound12", "Compound13", "Compound14", "Compound15"), 
  B2M = c(10,10,10,10,10,10,10,3,4,5,6,7,6,7,8),
HBD1 = c(10,10,10,10,10,10,10,10,10,10,6,7,6,7,8),
HBD2 = c(10,10,10,10,10,10,10,10,10,10,6,7,6,7,8),
HBD3 = c(10,10,10,10,10,10,10,10,10,10,6,7,6,7,8),
HBD4 = c(10,10,10,10,10,10,10,10,10,10,6,7,6,7,8),
LL37 = c(10,10,10,10,10,10,10,10,10,10,6,7,6,7,8),
REG3A = c(4,5,6,7,1,1,2,3,4,5,6,7,6,7,8),
PLA2G2A = c(4,5,6,7,5,1,2,3,4,5,6,7,6,7,8),
IL1B = c(1,3,2,2,2,2,1,10,4,5,2,2,2,2,2),
IL8 = c(1,3,2,2,2,2,1,10,4,5,2,2,2,2,2),
TNFalpha = c(1,3,2,2,2,2,1,10,4,5,2,2,2,2,2))

### tSNE 
# so i have 2 seperate tables diff_expression and compound_data
data_tsne <-tsne::tsne(compound_data[1:ncol(compound_data)-1])
data_tsne <-as.data.frame(data_tsne)
ggplot(data_tsne, aes(V1, V2)) +
  geom_point()

###
expression_tsne <-tsne::tsne(diff_expression[2:ncol(diff_expression)])
expression_tsne <- as.data.frame(expression_tsne)
ggplot(expression_tsne, aes(V1, V2, color = Compound)) +
  geom_point()

### An additional table with expression values + compound similarity score
# just merge the two tables! 
rownames(compound_data) <- rownames(diff_expression$Compound)
# compound_data$Compound
diff_expression
# diff_expression$Compound

#### Similarity Matrix
compound_DE_expression <-full_join(diff_expression,compound_data, by = "Compound")

merged_tsne <- tsne::tsne(compound_DE_expression[2:12])
ggplot(merged_tsne, aes(V1,V2)) ## okay, dropping, probably wont work

# transforming 64bit vector to fingerprint class
finger <-decodeCACTVS(lc50_sub3$Fingerprint2D[1])
finger2 <-decodeCACTVS(lc50_sub3$Fingerprint2D[2])

fingerprint::distance(finger, finger2)

# similarity matrix
fingerprints_list <- list()
for (i in 1:length(fingerprints[[1]])){
  print(fingerprints[[1]][i])
  fingerprints_list[i] <-decodeCACTVS(fingerprints[[1]][i])
}

similarity_matrix <-fp.sim.matrix(fingerprints_list)
rownames(similarity_matrix) <- diff_expression$Compound
colnames(similarity_matrix)  <- diff_expression$Compound
#frankly not sure how to go about it except for representing it as a heatmpap

transformed_matrix <-similarity_matrix %>% 
  as.data.frame() %>%
  rownames_to_column(.,"Compound") %>%
  pivot_longer(-c(Compound), names_to = "compound2", values_to = "similarity")  ##didnt check these functions completely, not even sure all of that is useful

similarity_expression <-full_join(transformed_matrix,diff_expression, by = "Compound" )
tsne_result <-as.data.frame(tsne::tsne(similarity_expression[4:ncol(similarity_expression)])) ##okay so that's actually not good, every value is duplicated, but if I could somehow remove them. Maybe instead do hierarchical cluster 

tsne_result$similarity <- similarity_expression$similarity

ggplot(tsne_result, aes(V1,V2, color = similarity)) +
  geom_point()

###### The final protocol

compound <-get_cid("Glyphosate", from = "name", match = "first", domain = "substance")

compound_info <- pc_prop(compound$cid, properties = c("IUPACName", "XLogP", "Fingerprint2D"))

binary2 <- list()
for (i in 1:nrow(compound_info)) {
  binary2[i] <- as.character(decodeCACTVS(compound_info$Fingerprint2D[i]))
}

matrix_glyphosate <-matrix(0,nrow = 15, ncol = 881) # or whatever length it'll be 
for (i in 1:length(binary)){
  matrix_glyphosate[i,] <- as.numeric(t(as.data.frame(strsplit(binary2[[i]], ""))))
}

#### Hierarchical cluster!!


