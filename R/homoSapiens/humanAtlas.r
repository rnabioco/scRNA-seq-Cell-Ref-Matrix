library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
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
humanGenes <- read_tsv(humanGenesFile)

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

common <- Reduce(intersect, list(rownames(GSE129933), rownames(GSE137710Melanoma), rownames(GSE137710Spleen), rownames(GSE147405)))
GSE129933[common,] # give you common rows in data frame 1  
GSE137710Spleen[common,] # give you common rows in data frame 2
GSE137710Melanoma[common,]
GSE147405[common,]
head(common)

#GSE129933 <- subset(GSE129933, rownames(GSE129933) %in% c(common))
#GSE137710Melanoma <- subset(GSE137710Melanoma, rownames(GSE137710Melanoma) %in% c(common))
#GSE137710Spleen <- subset(GSE137710Spleen, rownames(GSE137710Spleen) %in% c(common))
#GSE147405 <- subset(GSE147405, rownames(GSE147405) %in% c(common))

GSE129933 <- GSE129933[common, ]
GSE137710Melanoma <- GSE137710Melanoma[common, ]
GSE137710Spleen <- GSE137710Spleen[common, ]
GSE147405 <- GSE147405[common, ]
GSE147405 <- as.data.frame(GSE147405)



appendGenes <- function(humanGenes, newGSEFile)
{
  rownamesHumanGenes <- rownames(humanGenes)
  rownamesNewGSEFile <- rownames(newGSEFile)
  
  rowCountHumanGenes <- rowCount(humanGenes)
  rowCountNewGSEFile <- rowCount(newGSEFile)
  
  for (i in rowCountNewGSEFile)
  {
    
  }
}

humanAtlas <- bind_rows(GSE129933, GSE137710Melanoma, GSE137710Spleen, GSE147405, .id = NULL)
saveRDS(humanAtlas, "HumanAtlas.rds")
