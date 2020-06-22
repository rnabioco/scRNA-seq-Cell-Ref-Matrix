library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)

GSE113049Filename <- file.choose()
GSE113049 <- readRDS(GSE113049Filename)

GSE124952Filename <- file.choose()
GSE124952 <- readRDS(GSE124952Filename)

GSE143405Filename <- file.choose()
GSE143405 <- readRDS(GSE143405Filename)

GSE137710Filename <- file.choose()
GSE137710 <- readRDS(GSE137710Filename)