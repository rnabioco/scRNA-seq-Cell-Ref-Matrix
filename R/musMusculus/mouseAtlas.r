library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(readr)
library(digest)

GSE113049Filename <- file.choose()
GSE113049 <- readRDS(GSE113049Filename)

GSE124952Filename <- file.choose()
GSE124952 <- readRDS(GSE124952Filename)

GSE137710Filename <- file.choose()
GSE137710 <- readRDS(GSE137710Filename)

GSE143435D0Filename <- file.choose()
GSE143435_D0 <- readRDS(GSE143435D0Filename)

GSE143435D2Filename <- file.choose()
GSE143435_D2 <- readRDS(GSE143435D2Filename)

GSE143435D5Filename <- file.choose()
GSE143435_D5 <- readRDS(GSE143435D5Filename)

GSE143435D7Filename <- file.choose()
GSE143435_D7 <- readRDS(GSE143435D7Filename)

mouseGenesFile <- file.choose()
fullMouseGenes <- read_tsv(mouseGenesFile)

rm(GSE113049Filename)
rm(GSE124952Filename)
rm(GSE137710Filename)
rm(GSE143435D0Filename)
rm(GSE143435D2Filename)
rm(GSE143435D5Filename)
rm(GSE143435D7Filename)

common <- Reduce(intersect, list(rownames(GSE113049), rownames(GSE124952), rownames(GSE137710), rownames(GSE143435_D0), rownames(GSE143435_D2), rownames(GSE143435_D5), rownames(GSE143435_D7)))
GSE113049 <- GSE113049[common,]
GSE124952 <- GSE124952[common,]
GSE137710 <- GSE137710[common,]
GSE143435_D0 <- GSE143435_D0[common,]
GSE143435_D2 <- GSE143435_D2[common,]
GSE143435_D5 <- GSE143435_D5[common,]
GSE143435_D7 <- GSE143435_D7[common,]

GSE113049 <- as.data.frame(GSE113049)
GSE124952 <- as.data.frame(GSE124952)
GSE137710 <- as.data.frame(GSE137710)
GSE143435_D0 <- as.data.frame(GSE143435_D0)
GSE143435_D2 <- as.data.frame(GSE143435_D2)
GSE143435_D5 <- as.data.frame(GSE143435_D5)
GSE143435_D7 <- as.data.frame(GSE143435_D7)

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

GSE113049NewRefMatrix <- appendGenes(fullMouseGenes, GSE113049)

mouseAtlas <- bind_rows(GSE113049, GSE124952, GSE137710, GSE143435_D0, GSE143435_D2, GSE143435_D5, GSE143435_D7, .id = NULL)
saveRDS(mouseAtlas, "MouseAtlas.rds")
