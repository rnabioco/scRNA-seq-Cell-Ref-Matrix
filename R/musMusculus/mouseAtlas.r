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

GSE137710Filename <- file.choose()
GSE137710 <- readRDS(GSE137710Filename)

GSE143405D0Filename <- file.choose()
GSE143405_D0 <- readRDS(GSE143405D0Filename)

GSE143435D2Filename <- file.choose()
GSE143435_D2 <- readRDS(GSE143435D2Filename)

GSE143435D5Filename <- file.choose()
GSE143435_D5 <- readRDS(GSE143435D5Filename)

GSE143435D7Filename <- file.choose()
GSE143435_D7 <- readRDS(GSE143435D7Filename)