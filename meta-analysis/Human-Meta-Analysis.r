library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

humanAtlasFilename <- file.choose()
humanAtlas <- readRDS(humanAtlasFilename)

humanMetaAnalysis <- CreateSeuratObject(counts = humanAtlas, project = "Human-Meta-Analysis", min.cells = 3, min.features = 200)
humanMetaAnalysis
rm(humanAtlas)
gc()