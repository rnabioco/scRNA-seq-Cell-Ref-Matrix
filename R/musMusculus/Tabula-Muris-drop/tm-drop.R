library(Seurat)
library(clustifyr)
library(tidyverse)
library(here)

#get project root
proj_dir <- here()

tmp_dir <- file.path(proj_dir, "tmp_tm_files")
dir.create(tmp_dir, showWarnings = FALSE)

# warning, 14Gb download
download.file("https://ndownloader.figshare.com/articles/5821263/versions/3",
              destfile = file.path(tmp_dir, "tm.zip"))
unzip(file.path(tmp_dir, "tm.zip"), exdir = tmp_dir)

# find .Robj files (from 10x genomics drop experiments)
files_drop <- list.files(tmp_dir,
                         pattern = "droplet.*\\.Robj$",
                         full.names = TRUE)

# function to load .Robj file and assign to custom variable
load_obj <- function(f)
{
  env <- new.env()
  nm <- load(f, env)[1]
  env[[nm]]
}

for (n in 1:length(files_drop)) {
  print(files_drop[n])
  # load and update seurat object to v3
  s <- load_obj(files_drop[n]) %>%
    UpdateSeuratObject(.)

  # average two levels of cell type classifications
  avg <- object_ref(s, cluster_col = "cell_ontology_class")
  avg2 <- object_ref(s, cluster_col = "free_annotation")
  avg <- cbind(avg, avg2)

  # add tissue to end of column name
  colnames(avg) <- str_c(colnames(avg),
                         str_sub(str_extract(basename(files_drop[n]), "_.+?_"),
                                 2,
                                 -2),
                         sep = "-")
  if (n==1) {
    ref <- avg
  } else {
    ref <- cbind(ref, avg)
  }
}

# put output into ref_matrices directory
out_dir <- file.path(proj_dir, "ref_matrices", "musMusculus", "Tabula-Muris-drop")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(ref, file.path(out_dir, "Tabula-Muris-drop.rds"))

# delete tmp dir with large R objects
#unlink(tmp_dir, recursive = TRUE, force = TRUE)
