# conda env: r-mofa-singler 
# this is for use part of the data as reference and do the supervised learning:

library(Seurat)
library(SeuratDisk)
library(SingleR)
       
train_input_dir = "/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TrainData/"
output_dir = '/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TestData/'

train_singler = function(train_file="seu_obj_few25_s400.h5ad", test_file='seu_obj_ForTest.h5ad') {
       train_pre = strsplit(train_file, split='.',fixed=T)[[1]][1]
       test_pre = strsplit(test_file, split='.',fixed=T)[[1]][1]
       
       Convert(paste0(train_input_dir, train_file), dest = "h5seurat", overwrite = TRUE)
       s400_ob <- LoadH5Seurat(paste0(train_input_dir, train_pre, ".h5seurat"), meta.data = FALSE, misc = FALSE)
       
       s400_ob_mat = s400_ob@assays$logcounts@data
       
       #Convert(paste0(input_dir, "TestData/","seu_obj_ForTest.h5ad"), dest = "h5seurat", overwrite = TRUE)
       #test_ob <- LoadH5Seurat(paste0(input_dir, "TestData/","seu_obj_ForTest.h5ad"), meta.data = FALSE, misc = FALSE, commands= FALSE,tools = FALSE )
       
       # load test data from csv:
       test_mat = read.csv(paste0(output_dir, test_pre, "_df.csv"), row.names = 1, header= TRUE)
       test_mat_m = as.matrix(test_mat)
       
       # load label for the s400_ob_mat:
       s400_ob_mat_label = read.csv(paste0(train_input_dir,train_pre,".csv"))
       
       # run singleR and output the labels
       
       pred.test <- SingleR(test=test_mat_m, ref=s400_ob_mat, labels=s400_ob_mat_label$cell_type, de.method="wilcox", de.n=50)
       table(pred.test$labels)
       write.table(pred.test['labels'], file=paste(output_dir ,test_pre,"_", train_pre, "_singler_labels",".csv",sep=""), row.names = TRUE, col.names = NA, sep=",", quote = FALSE)
       
}

#train_singler(train_file="seu_obj_few25_s400.h5ad", test_file='seu_obj_ForTest.h5ad')

for (train_file in list.files(path = train_input_dir)) {
       
       if (strsplit(train_file, split='.',fixed=T)[[1]][2] == "h5ad") {
              print(train_file)
              train_singler(train_file,test_file='seu_obj_ForTest.h5ad')
       }
}



# output the test df from h5ad, since above loading has issue.
#import scanpy as sc
#test = sc.read_h5ad("/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TestData/seu_obj_ForTest.h5ad")
#import pandas as pd 
#test_df = pd.DataFrame(data=test.layers['logcounts'].toarray(), index=test.obs_names, columns=test.var_names)
#test_df_t = test_df.T
#test_df_t.to_csv("/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TestData/test_df.csv")
