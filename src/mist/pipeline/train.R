#-------------------------------------------------------------------------------
# PipeElement: Extreme Gradient Boosting Algorithm -> Training
#-------------------------------------------------------------------------------
# --- AUX
# infolder <- "./outputs/r-objects/"
# epochs           <- 150
# e                <- 0.3
# g                <- 0
# md               <- 6
# mcwt             <- 1
# ss               <- 1
# csbt             <- 1
# v <- 2

Sys.setenv(CUDA="11.1")

#------------------ Requirements, dependencies and libraries
pkgs <- c('pacman','dplyr','docopt','xgboost', 'data.table', 'ggplot2', 'caret')        # Packages list
gitpkgs <- "agnesdeng/mixgb"                                                           # GitHub packages list

suppressMessages(if (!require("BiocManager", character.only = TRUE)) {                           #  
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

#---  EXPERIMENTAL --- #
# ADD REPO: https://launchpad.net/~c2d4u.team/+archive/ubuntu/c2d4u4.0+
# Install dependencies: RcppGSL and RcppZiggurat
# info  <- devtools::package_info(gitpkgs) # Verify if github packages are installed and/or handle their installation
# if(is.na(info$ondiskversion)){
#   paste0(gitpkgs, " package is missing. Trying to install from github...")
#   devtools::install_github("agnesdeng/mixgb")
# }else {
#   print("Github packages are installed and ready to use!")
# }
#--- /EXPERIMENTAL --- #

"Mesothelial Cancer Cell Line Data Analysis (M.C.Li.D.A) Pipeline - Training

Usage: CLI_MiST-XGBoost.R --infolder=<folder> --epochs=<value> --max_depth=<value> --eta=<value> --gamma=<value> --colsample_bytree=<value> --min_child_weight=<value> --subsample=<value> --verbose=<value>

Options:
  -h --help                  Show this screen.
  --infolder=<folder>        Folder where the outputs are placed.
  --epochs=<value>           Number of training epochs for the model.
  --max_depth=<value>
  --eta=<value>
  --gamma=<value>
  --colsample_bytree=<value>
  --min_child_weight=<value>
  --subsample=<value>
  --verbose=<value>          If set to 1 prints all messages.
  --version                  
"-> doc

print("Loading required packages...")
library(pacman)
pacman::p_load(pkgs, install = TRUE, character.only = TRUE)
pacman::p_loaded()

arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

#------------------Functions
source("./pipeline/functions.R")

#------------------ Load dataset and parameters into R environment 
infolder         <- arguments$infolder
epochs           <- as.numeric(arguments$epochs)
md               <- as.numeric(arguments$max_depth)
e                <- as.numeric(arguments$eta)           
g                <- as.numeric(arguments$gamma)
csbt             <- as.numeric(arguments$colsample_bytree)
mcwt             <- as.numeric(arguments$min_child_weight)
ss               <- as.numeric(arguments$subsample)   
verbose                <- as.numeric(arguments$verbose)

print(paste0("Running XGBoost algorithm with inputs and labels from  ", infolder))

#------------------ File paths to maintain the file system integrity
mname   <- "MiST-XGBoost.model"
paths <- list(pdf       = "./outputs/pdf/", 
              models    = "./outputs/models/", 
              m.eval    = "./outputs/models/model_eval/", 
              r.objects = "./outputs/r-objects/"
              )
fnames <- list(model = "MiST-XGBoost.model",
               DTrain = "train.Dmatrix",
               DTest  = "test.Dmatrix",
               learning.set = "learning.set"
               )

files <- list(
  ifile  = paste0(infolder, "inputs.RDS"),
  lfile  = paste0(infolder, "labels.RDS"),
  idfile = paste0(infolder, "sid.RDS"),
  dfile  = paste0(infolder, "dataframe.RDS")
)

#------------------ Read the saved R objects
pipe.data <- list(
  inputs    = readRDS(files$ifile),
  labels    = readRDS(files$lfile),
  sid       = readRDS(files$idfile),
  dataframe = readRDS(files$dfile)
)

#------------------ Sample the dataset (70%-30%)
set.seed(312)
qty <- 0.7
learning.set <- sampleDataset(qty,pipe.data$dataframe)

#----------------- Impute missing values (not working yet)
# heatmap(as.matrix(learning.set$train[,-c(1:2)]), Colv = learning.set$train$labels)

#------------------ Create XGB matrices
train.matrix <- xgb.DMatrix(data = as.matrix(learning.set$train[,-c(1:2)]), 
                            label = learning.set$train$labels, missing = 99) # Convert into a dense matrix, set NA's as missing values to be imputed by xgboost
test.matrix  <- xgb.DMatrix(data = as.matrix(learning.set$test[,-c(1:2)]), 
                            label = learning.set$test$labels, missing = 99)

#------------------ Set the parameters
params    <- list(booster         = "gbtree",
                 eval_metric      = "mlogloss",
                 objective        = "multi:softmax",
                 num_class        = 5,
                 eta              = e,
                 gamma            = g,
                 min_child_weight = mcwt,
                 max_depth        = md,
                 colsample_bytree = csbt,
                 subsample        = ss)

save_name <- "XGBoost_logloss.model"
seed      <- 314
watchlist <- list(train = train.matrix, test = test.matrix)

#------------------ Run Extreme Gradient Boosting algorithm (XGBOOST)
MiST_model <- createModel(data = train.matrix, 
                          params = params, 
                          watchlist = watchlist, 
                          epochs = epochs, 
                          save_name = save_name, 
                          verbose = verbose, 
                          seed = seed)
delta <-list( d.basic = (MiST_model$evaluation_log[tail(MiST_model$evaluation_log$iter, 1),]$test_mlogloss - 
                           MiST_model$evaluation_log[tail(MiST_model$evaluation_log$iter, 1),]$train_mlogloss)) # 

e.log <- data.frame(MiST_model$evaluation_log)                # Extract the xgboost training evaluation log
imatrix <- xgb.importance(colnames(learning.set$train[,-c(1:2)]), MiST_model) # Extract the matrix of feature importance


#------------------ Evaluation plot
pdf(paste0(paths[1],"summary.pdf"), height = 11.69, width = 8.27)
plots <- p(e.log, MiST_model, paths, as.matrix(learning.set$train[,-c(1:2)]), imatrix)
dev.off()

jpeg(paste0(paths[1],"multi_trees.pdf"))
xgb.plot.multi.trees(MiST_model, feature_names = colnames(learning.set$test[,-c(1:2)]),features_keep = 5,render = TRUE)
dev.off()

#------------------ Cross-validate the model
set.seed(314)
cv <- xgb.cv(params = params, data = train.matrix, nfold = 5, showsd = T, stratified = T, print_every_n = 5, early_stopping_rounds = 20, nrounds = 100, verbose = verbose)
lv <- crossVal(cv)
pCV(paths[3], cv$evaluation_log$test_mlogloss_mean)
#------------------ Tune the parameters
epochs <- lv
params <- list(booster          = "gbtree",
               eval_metric      = "mlogloss", 
               objective        = "multi:softmax", 
               num_class        = 5, 
               eta              = e, 
               gamma            = g, 
               min_child_weight = mcwt, 
               max_depth        = md, 
               colsample_bytree = csbt, 
               subsample        = ss)

#------------------ Re-train the model using updated parameters
save_name      <- "XGBoost_logloss.model"
MiST_model.re  <- createModel(data = train.matrix, 
                              params = params, 
                              watchlist = watchlist, 
                              epochs = epochs, 
                              save_name = save_name,
                              verbose = verbose,
                              seed = seed)

delta$d.re     <- (MiST_model.re$evaluation_log[tail(MiST_model.re$evaluation_log$iter, 1),]$test_mlogloss - 
                     MiST_model.re$evaluation_log[tail(MiST_model.re$evaluation_log$iter, 1),]$train_mlogloss) # Add the new change in log-loss

if(delta$d.re < delta$d.basic){
  paste0("Successefully reduced train-test error difference by -", delta$d.basic - delta$d.re)
}



#------------------ Evaluation plot
e.cv.re <- data.frame(MiST_model.re$evaluation_log)

pCV(paths[3],e.cv.re = e.cv.re)

# pdf(paste0(paths[1],"mlogloss.pdf")) # Save the plot as pdf
# scatter.smooth(e.re, xlab = "No. of Rounds", ylab = "Change in mLogLoss")
# plot(e.re$iter, e.re$train_mlogloss, col = 'blue')
# lines(e.re$iter, e.re$test_mlogloss, col = 'red')
# dev.off()

#------------------ Save generated outputs
model <- xgb.save(MiST_model.re, paste0(paths[2],fnames[1]))
xgb.DMatrix.save(train.matrix, paste0(paths[4],fnames[2]))
xgb.DMatrix.save(test.matrix, paste0(paths[4],fnames[3]))
saveRDS(learning.set, paste0(paths[4],fnames[4]))

if(file.exists(paste0(paths[2],fnames[1]))){
  print(paste0("Model were succesefully saved as '", fnames[1], "'."))
}
paste0("Training and testing matrices are saved as xgb.Dmatrix files.")
paste0("Learning set was saved as learning.set.")

