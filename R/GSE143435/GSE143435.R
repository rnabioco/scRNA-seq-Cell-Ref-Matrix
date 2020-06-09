library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

D0_FACSatlas <- read_csv("~/Reference-Matrix-Generation/data/GSE143435/GSE143435.csv")
D0_FACSatlas <- D0_FACSatlas %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  as.matrix() %>%
  t()
D0_FACSatlas[1:5, 1:5]

D0_FACSatlasMetadata <- read_text("~/Reference-Matrix-Generation/data/GSE143435/GSE143435.csv")
sum(colnames(D0_FACSatlas) %in% D0_FACSatlasMetadata$)
ncol(D0_FACSatlas)