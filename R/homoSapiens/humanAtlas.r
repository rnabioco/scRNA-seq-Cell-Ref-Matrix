library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(readr)
library(digest)

GSE129933Filename <- file.choose()
GSE129933Matrix <- readRDS(GSE129933Filename)

GSE137710FilenameMelanoma <- file.choose()
GSE137710Melanoma <- readRDS(GSE137710FilenameMelanoma)

GSE137710FilenameSpleen <- file.choose()
GSE137710Spleen <- readRDS(GSE137710FilenameSpleen)

GSE147405Filename <- file.choose()
GSE147405Matrix <- readRDS(GSE147405Filename)

humanGenesFile <- file.choose()
humanGenesTSV <- read_tsv(humanGenesFile)
fullHumanGenes <- as.data.frame(humanGenesTSV)
rm(humanGenesTSV)
humanGenesVector <- as.vector(fullHumanGenes[,1])

GSE129933 <- as.data.frame(GSE129933Matrix)
GSE147405 <- as.data.frame(GSE147405Matrix)
GSE137710Melanoma <- as.data.frame(GSE137710Melanoma)
GSE137710Spleen <- as.data.frame(GSE137710Spleen)

rm(GSE129933Matrix)
rm(GSE147405Matrix)
rm(GSE129933Filename)
rm(GSE147405Filename)
rm(GSE137710FilenameMelanoma)
rm(GSE137710FilenameSpleen)
rm(humanGenesFile)

appendGenes <- function(humanGenesVector, GSEMatrix)
{
  rownamesGSEMatrix <- rownames(GSEMatrix) #Get rownames from GSEMatrix (new GSE file)
  
  rowCountHumanGenes <- nrow(humanGenesVector) #Calculate number of rows from list of full human genes
  rowCountNewGSEFile <- nrow(GSEMatrix) #Calculate number of rows of GSE matrix
  
  missing_rows <- setdiff(humanGenesVector, rownamesGSEMatrix) #Use setdiff function to figure out rows which are different/missing from GSE matrix
  missing_rows #Display missing rows
  
  zeroExpressionMatrix <- matrix(0, nrow = length(missing_rows), ncol = ncol(GSEMatrix)) #Create a placeholder matrix with zeroes and missing_rows length as row length
  #zeroExpressionMatrix[1:3, 1:3] #Check first three entries of zeroExpressionMatrix
  dim(zeroExpressionMatrix) #Check dimensions of matrix
  
  rownames(zeroExpressionMatrix) <- missing_rows #Assign row names
  colnames(zeroExpressionMatrix) <- colnames(GSEMatrix) #Assign column names
  #zeroExpressionMatrix[1:3, 1:3] #Check col and row names assigned correctly
  
  fullMatrix <- rbind(GSEMatrix, zeroExpressionMatrix) #Bind GSEMatrix and zeroExpressionMatrix together
  dim(fullMatrix) #Check dimensions of fullMatrix
  
  #Reorder matrix
  fullMatrix <- fullMatrix[humanGenesVector, ] #Reorder fullMatrix to preserve gene order
  #fullMatrix[1:10, 1:3] #Check reordering
  return(fullMatrix) #Return fullMatrix
}

GSE129933Matrix <- as.matrix(GSE129933)
GSE129933NewRefMatrix <- appendGenes(humanGenesVector = humanGenesVector, GSEMatrix = GSE129933Matrix)
head(GSE129933NewRefMatrix)

GSE137710MelanomaMatrix <- as.matrix(GSE137710Melanoma)
GSE137710NewMelanomaRefMatrix <- appendGenes(humanGenesVector = humanGenesVector, GSEMatrix = GSE137710MelanomaMatrix)
head(GSE137710NewMelanomaRefMatrix)

GSE137710SpleenMatrix <- as.matrix(GSE137710Spleen)
GSE137710NewSpleenRefMatrix <- appendGenes(humanGenesVector = humanGenesVector, GSEMatrix = GSE137710SpleenMatrix)
head(GSE137710NewSpleenRefMatrix)

GSE147405Matrix <- as.matrix(GSE147405)
GSE147405NewRefMatrix <- appendGenes(humanGenesVector = humanGenesVector, GSEMatrix = GSE147405Matrix)
head(GSE147405NewRefMatrix)

humanAtlas <- cbind(GSE129933NewRefMatrix, GSE137710NewMelanomaRefMatrix, GSE137710NewSpleenRefMatrix, GSE147405NewRefMatrix, .id = NULL)

colnames(humanAtlas) <- c("Erythrocytes (GSE129933)",
                          "Fibroblasts (GSE129933)",
                          "Hepatocytes (GSE129933)",
                          "LEC (GSE129933)",
                          "Lymphocytes (GSE129933)",
                          "Monocytes (GSE129933)",
                          "PEC (GSE129933)",
                          "B cells (GSE137710 Melanoma)",
                          "melanoma (GSE137710 Melanoma)",
                          "Myeloid (GSE137710 Melanoma)",
                          "T/NK (GSE137710 Melanoma)",
                          "AS DC (GSE137710 Spleen)",
                          "CCR7+ cDC2 (GSE137710 Spleen)",
                          "cDC1 (GSE137710 Spleen)",
                          "CLEC10A- cDC2 (GSE137710 Spleen)",
                          "CLEC10A+ cDC2 (GSE137710 Spleen)",
                          "Mitotic cDC1 (GSE137710 Spleen)",
                          "Mitotic cDC2 (GSE137710 Spleen)",
                          "OVCA420 (GSE147405)"
                          )
saveRDS(humanAtlas, "HumanAtlas.rds")
