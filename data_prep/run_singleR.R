# this is for use part of the data as reference and do the supervised learning:

library(Seurat)
library(SeuratDisk)
       
input_dir = "/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/"


Convert(paste0(input_dir, "TrainData/","seu_obj_few25_s400.h5ad"), dest = "h5seurat", overwrite = TRUE)
s400_ob <- LoadH5Seurat(paste0(input_dir, "TrainData/","seu_obj_few25_s400.h5seurat"))


