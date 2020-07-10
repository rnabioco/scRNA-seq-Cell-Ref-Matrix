library(dplyr)
library(Seurat)
library(clustifyr)
library(tidyverse)
library(readr)
library(digest)
library(here)

# figure out project root
proj_dir <- here()

# get scripts
source(file.path(proj_dir, "R", "utils", "utils.r"))

# path to matrices
ref_matrix_dir <- file.path(proj_dir, "ref_matrices", "musMusculus")

# find all .rds files in the ref_matrices directory,
# name vector with filename without .rds
ref_matrices_fns <- list.files(ref_matrix_dir,
                               pattern = ".rds$",
                               full.names = TRUE,
                               recursive = TRUE)
names(ref_matrices_fns) <- basename(ref_matrices_fns) %>% str_remove(".rds$")

# remove studies in records_to_drop
records_to_drop <- c("GSE137710MouseSpleen")
idx_to_keep <- which(!names(ref_matrices_fns) %in% records_to_drop)
ref_matrices_fns <- ref_matrices_fns[idx_to_keep]

# path to mouse genes
mouse_genes_fn <- file.path(proj_dir,
                            "data",
                            "geneList",
                            "mouse_genes.tsv.gz")

# output filename for atlas
atlas_dir <- file.path(proj_dir, "atlas", "musMusculus")
dir.create(atlas_dir, showWarnings = FALSE, recursive = TRUE)
atlas_fn <- file.path(atlas_dir, "MouseAtlas.rds")

build_atlas(ref_matrices_fns, mouse_genes_fn, atlas_fn)
