library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

mat_FibroblastNucleiDiabetes <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148946/suppl/GSE148946_nuclei_normdata.csv.gz")
mat_FibroblastNucleiDiabetes <- mat_FibroblastNucleiDiabetes %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  as.matrix() %>%
  t()
mat_FibroblastNucleiDiabetes[1:5, 1:5]

meta_FibroblastNucleiDiabetes <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148946/suppl/GSE148946_nuclei_metadata.csv.gz")
sum(colnames(mat_FibroblastNucleiDiabetes) %in% meta_FibroblastNucleiDiabetes$celltype)
ncol(mat_FibroblastNucleiDiabetes)

FibroblastNucleiDiabetes <- CreateSeuratObject(counts = mat_FibroblastNucleiDiabetes, project = "FibroblastNucleiDiabetes", min.cells = 3, min.features = 200)
FibroblastNucleiDiabetes
