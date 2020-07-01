library(dplyr)
library(Seurat)
library(clustifyr)
library(tidyverse)
library(readr)
library(digest)

GSE113049 <- readRDS(file.path("~/Reference-Matrix-Generation/ref_matrices/musMusculus/GSE113049/GSE113049.rds"))

GSE124952 <- readRDS(file.path("~/Reference-Matrix-Generation/ref_matrices/musMusculus/GSE124952/GSE124952.rds"))

GSE137710 <- readRDS(file.path("~/Reference-Matrix-Generation/ref_matrices/musMusculus/GSE137710MouseCopy/GSE137710MouseSpleen.rds"))

GSE143435_D0 <- readRDS(file.path("~/Reference-Matrix-Generation/ref_matrices/musMusculus/GSE143435/GSE143435D0.rds"))

GSE143435_D2 <- readRDS(file.path("~/Reference-Matrix-Generation/ref_matrices/musMusculus/GSE143435/GSE143435D2.rds"))

GSE143435_D5 <- readRDS(file.path("~/Reference-Matrix-Generation/ref_matrices/musMusculus/GSE143435/GSE143435D5.rds"))

GSE143435_D7 <- readRDS(file.path("~/Reference-Matrix-Generation/ref_matrices/musMusculus/GSE143435/GSE143435D7.rds"))

mouseGenesFile <- file.choose()
mouseTSV <- read_tsv(mouseGenesFile)
fullMouseGenes <- as.data.frame(mouseTSV)
rm(mouseTSV)
mouseGenesVector <- as.vector(fullMouseGenes[,1])

rm(GSE113049Filename)
rm(GSE124952Filename)
rm(GSE137710Filename)
rm(GSE143435D0Filename)
rm(GSE143435D2Filename)
rm(GSE143435D5Filename)
rm(GSE143435D7Filename)

appendGenes <- function(mouseGenesVector, GSEMatrix)
{
  rownamesGSEMatrix <- rownames(GSEMatrix) #Get rownames from GSEMatrix (new GSE file)
  
  rowCountHumanGenes <- nrow(mouseGenesVector) #Calculate number of rows from list of full human genes
  rowCountNewGSEFile <- nrow(GSEMatrix) #Calculate number of rows of GSE matrix
  
  missing_rows <- setdiff(mouseGenesVector, rownamesGSEMatrix) #Use setdiff function to figure out rows which are different/missing from GSE matrix
  missing_rows #Display missing rows
  
  zeroExpressionMatrix <- matrix(0, nrow = length(missing_rows), ncol = ncol(GSEMatrix)) #Create a placeholder matrix with zeroes and missing_rows length as row length
  zeroExpressionMatrix[1:3, 1:3] #Check first three entries of zeroExpressionMatrix
  dim(zeroExpressionMatrix) #Check dimensions of matrix
  
  rownames(zeroExpressionMatrix) <- missing_rows #Assign row names
  colnames(zeroExpressionMatrix) <- colnames(GSEMatrix) #Assign column names
  zeroExpressionMatrix[1:3, 1:3] #Check col and row names assigned correctly
  
  fullMatrix <- rbind(GSEMatrix, zeroExpressionMatrix) #Bind GSEMatrix and zeroExpressionMatrix together
  dim(fullMatrix) #Check dimensions of fullMatrix
  
  #Reorder matrix
  fullMatrix <- fullMatrix[mouseGenesVector, ] #Reorder fullMatrix to preserve gene order
  fullMatrix[1:10, 1:3] #Check reordering
  return(fullMatrix) #Return fullMatrix
}

GSE113049Matrix <- as.matrix(GSE113049)
GSE113049NewRefMatrix <- appendGenes(mouseGenesVector = mouseGenesVector, GSEMatrix = GSE113049Matrix)
head(GSE113049NewRefMatrix)

GSE124952Matrix <- as.matrix(GSE124952)
GSE124952NewRefMatrix <- appendGenes(mouseGenesVector = mouseGenesVector, GSEMatrix = GSE124952Matrix)
head(GSE124952NewRefMatrix)

GSE137710Matrix <- as.matrix(GSE137710)
GSE137710NewRefMatrix <- appendGenes(mouseGenesVector = mouseGenesVector, GSEMatrix = GSE137710Matrix)
head(GSE137710NewRefMatrix)

GSE143435_D0Matrix <- as.matrix(GSE143435_D0)
GSE143435_D0NewRefMatrix <- appendGenes(mouseGenesVector = mouseGenesVector, GSEMatrix = GSE143435_D0Matrix)
head(GSE143435_D0NewRefMatrix)

GSE143435_D2Matrix <- as.matrix(GSE143435_D2)
GSE143435_D2NewRefMatrix <- appendGenes(mouseGenesVector = mouseGenesVector, GSEMatrix = GSE143435_D2Matrix)
head(GSE143435_D2NewRefMatrix)

GSE143435_D5Matrix <- as.matrix(GSE143435_D5)
GSE143435_D5NewRefMatrix <- appendGenes(mouseGenesVector = mouseGenesVector, GSEMatrix = GSE143435_D5Matrix)
head(GSE143435_D5NewRefMatrix)

GSE143435_D7Matrix <- as.matrix(GSE143435_D7)
GSE143435_D7NewRefMatrix <- appendGenes(mouseGenesVector = mouseGenesVector, GSEMatrix = GSE143435_D7Matrix)
head(GSE143435_D7NewRefMatrix)

mouseAtlas <- cbind(GSE113049NewRefMatrix, GSE124952NewRefMatrix, GSE137710NewRefMatrix, GSE143435_D0NewRefMatrix, GSE143435_D2NewRefMatrix, GSE143435_D5NewRefMatrix, GSE143435_D7NewRefMatrix, .id = NULL)
saveRDS(mouseAtlas, "MouseAtlas.rds")

colnames(mouseAtlas) <- c("Basal (GSE113049)", 
                          "Ciliated = Ciliated (GSE113049)", 
                          "Club = Club (GSE113049)", 
                          "Endothelial/Fibroblast = Endothelial/Fibroblast (GSE113049)", 
                          "Injured AEC2: Cell Cycle Arrest (GSE113049)", 
                          "Injured AEC2: Proliferating (GSE113049)", 
                          "Injured AEC2: Transdifferentiating (GSE113049)", 
                          "Macrophage (GSE113049)", 
                          "Naive AEC1 (GSE113049)", 
                          "Naive AEC2 (GSE113049)", 
                          "Other Injured AEC2 (GSE113049)", 
                          "Astro = Astro (GSE124952)", 
                          "Endo = Endothelial (GSE124952)", 
                          "Excitatory = Excitatory (GSE124952)",
                          "Inhibitory = Inhibitory (GSE124952)",
                          "Microglia = Microglia (GSE124952)",
                          "NF Oligo (GSE124952)",
                          "Oligo (GSE124952)",
                          "OPC (GSE124952)",
                          "CCR7hiDC (GSE137710)",
                          "cDC1 (GSE137710)",
                          "cDC2 Mixed (GSE137710)",
                          "cDC Tbet- (GSE137710)",
                          "cDC Tbet+ (GSE137710)",
                          "Monocyte (GSE137710)", 
                          "Singlec - H DC (GSE137710)", 
                          "Endothelial (GSE143435_D0)",
                          "FAPs (GSE143435_D0)",
                          "Immune (GSE143435_D0)",
                          "Neural/Glial/Schwann cells (GSE143435_D0)",
                          "Platelets (GSE143435_D0)",
                          "Smooth muscle cells (GSE143435_D0)",
                          "Tenocytes (GSE143435_D0)",
                          "B cells (GSE143435_D0)",
                          "Endothelial (GSE143435_D2)",
                          "FAPs (GSE143435_D2)",
                          "Immune (GSE143435_D2)",
                          "Macrophages/APCs (GSE143435_D2)",
                          "MuSCs and progenitors (GSE143435_D2)",
                          "Anti-inflamatory macrophages (GSE143435_D5)",
                          "B cells (GSE143435_D5)",
                          "Endothelial (GSE143435_D5)",
                          "FAPs (GSE143435_D5)",
                          "Macrophages/APCs/T-cells (GSE143435_D5)",
                          "Mature skeletal muscles (GSE143435_D5)",
                          "MuSCs and progenitors (GSE143435_D5)",
                          "NK cells (GSE143435_D5)",
                          "Anti-inflamatory macrophages (GSE143435_D7)",
                          "Endothelial (GSE143435_D7)",
                          "FAPs (GSE143435_D7)",
                          "FAPs activated (GSE143435_D7)",
                          "Macrophages/APCs (GSE143435_D7)",
                          "Macrophages/APCs 2 (GSE143435_D7)", 
                          "Mature skeletal muscles (GSE143435_D7)",
                          "MuSCs and progenitors (GSE143435_D7)",
                          "NK/T cells (GSE143435_D7)",
                          "Tenocytes (GSE143435_D7)"
                          )