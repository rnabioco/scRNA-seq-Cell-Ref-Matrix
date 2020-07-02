library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

humanAtlas <- readRDS(file.path("~/Reference-Matrix-Generation/atlas/homoSapiens/HumanAtlas.rds"))

humanMetaAnalysis <- CreateSeuratObject(counts = humanAtlas, project = "Human-Meta-Analysis", min.cells = 3, min.features = 200)
humanMetaAnalysis
rm(humanAtlas)
gc()
