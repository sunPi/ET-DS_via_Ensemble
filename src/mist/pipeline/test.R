#-------------------------------------------------------------------------------
# PipeElement: Extreme Gradient Boosting Algorithm -> Predicting (Testing)
#-------------------------------------------------------------------------------
#------------------ Requirements, dependencies and libraries
# caTools
#--- AUX
# infolder <- "./outputs/"

pkgs <- c('pacman','docopt','dplyr','xgboost', 'pacman', 'data.table', 'writexl', 'caret')

suppressMessages(if (!require("BiocManager", character.only = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
} else {
  ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
  if (any(!ipkgs)) {
    BiocManager::install(pkgs[!ipkgs])
    install.packages(pkgs[!ipkgs])
  } else {
    message("\n\nCool! your machine has everything is needed.\n\n")
  }
})

print("Loading required packages...")
library(pacman)
pacman::p_load(pkgs, install = TRUE, character.only = TRUE)
pacman::p_loaded()

"Mesothelial Cancer Cell Line Data Analysis (M.C.Li.D.A) Pipeline - Prediction

Usage: CLI_MiST-XGBoost_predict.R --infolder=<folder> 

Options:
  -h --help                  Show this screen.
  --infolder=<folder>        Folder where the outputs are placed.
  --version                  Pipeline Pre-release version
"-> doc

arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

#------------------Functions
source("./pipeline/functions.R")

#------------------ Load dataset and parameters into R environment 
infolder <- arguments$infolder
paths    <- list(pdf       = "./outputs/pdf/", 
                 models    = "./outputs/models/", 
                 m.eval    = "./outputs/models/model_eval/", 
                 r.objects = "./outputs/r-objects/")

fnames   <- list(model        = "MiST-XGBoost.model",
                 sid          = "sid.RDS",
                 learning.set = "learning.set")

#------------------ Function that tests difference in gene features between training and test data frames
# gene_selection <- colnames(pi)
# str(test)
# tgenes <- colnames(test)
# m <- match(gene_selection, tgenes)
# ind <- which(tgenes %in% gene_selection)
# 
# true_genes_list <- c(tgenes[ind])
# true_genes_list

#------------------ Load the trained model
MiST.model   <- xgb.load(paste0(paths[2],fnames[1]))
sampleID     <- readRDS(paste0(paths[4],fnames[2]))
learning.set <- readRDS(paste0(paths[4],fnames[3]))
#------------------ Test the model
print(paste0("Predicting on patients via XGBoost with inputs from ", infolder))
MiST.predict <- predict(MiST.model, newdata = as.matrix(learning.set$test[,-c(1:2)]), missing = 99)
MiST.predict
learning.set$test$labels
# confusion.matrix <- confusionMatrix(as.factor(MiST.predict), as.factor(learning.set$test$labels))
metrics <- calcMetrics(MiST.predict, learning.set$test$labels)
metrics

#------------------ Save predictions
out <- data.frame(cbind(sampleID    = learning.set$test$sample, 
             response    = learning.set$test$labels, 
             predictions = MiST.predict
             ))

out[1,4] <- metrics$accuracy 
out[1,5] <- metrics$precision
out[1,6] <- metrics$recall
out[1,7] <- metrics$F1
colnames(out)[4:7] <- c("accuracy", "precision", "recall", "F1")

write.csv(out, paste0(infolder,"output.csv"))


