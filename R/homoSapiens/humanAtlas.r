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

GSE129933 <- as.data.frame(GSE129933Matrix)
GSE147405 <- as.data.frame(GSE147405Matrix)

rm(GSE129933Matrix)
rm(GSE147405Matrix)

bind_rows(GSE129933, GSE137710Melanoma, GSE137710Spleen, GSE147405, .id = NULL)