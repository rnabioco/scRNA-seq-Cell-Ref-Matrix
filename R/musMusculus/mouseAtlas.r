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

appendGenes <- function(mouseGenes, newGSEFile)
{
  rownamesNewGSEFile <- rownames(newGSEFile)
  
  rowCountHumanGenes <- nrow(mouseGenes)
  rowCountNewGSEFile <- nrow(newGSEFile)
  
  for (i in rowCountNewGSEFile)
  {
    searchGene <- rownamesNewGSEFile[i]
    for (j in rowCountHumanGenes)
    {
      if (searchGene != mouseGenes[j,1])
      {
        newGeneArray <- searchGene
      }
    }
  }
  return(newGeneArray)
}

appendGenes(fullMouseGenes, GSE113049)

mouseAtlas <- bind_rows(GSE113049, GSE124952, GSE137710, GSE143435_D0, GSE143435_D2, GSE143435_D5, GSE143435_D7, .id = NULL)
saveRDS(mouseAtlas, "MouseAtlas.rds")
