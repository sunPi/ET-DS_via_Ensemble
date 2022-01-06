#-------------------------------------------------------------------------------
# Evolutionary Trajectory-Drug Sensitivity via ENSEMBLE Learning Pipeline
#-------------------------------------------------------------------------------
# NUMERICAL CODING IN CANCER BIOSTATISTICAL ANALYSIS 
# 0   lack		
# -1	loss -> WE expect this to be  Associated with sensitivity. 	Bridge analogy for  Synthetic lethality as a response. 
# -2	deletion		
# 1	  gain		
# 2	  amplifications		
# Missing (99)	  upd	-> Unparental  Disomy  (mimics deletion/amplification)	
#
# Complete Response (CR): Disappearance of all target lesions
# Partial Response (PR): At least a 30% decrease in the sum of the LD of target lesions, taking as reference the baseline sum LD
# Progressive Disease (PD): At least a 20% increase in the sum of the LD of target lesions, taking as reference the smallest sum LD recorded since the treatment started or the appearance of one or more new lesions
# Stable Disease (SD): Neither sufficient shrinkage to qualify for PR nor sufficient increase to qualify for PD, taking as reference the smallest sum LD since the treatment started
#
# C  - clonal loss 2
# Cd - clonal_deletion -2
# Ca - clonal amplificaiton
# S  - subclonal loss 1
# Sd - subclonal deletion -1
# Sa - subclonal amplification
#
# XGBOOST Parameters
# - verbosity: Verbosity of printing messages. Valid values of 0 (silent), 1 (warning), 2 (info), and 3 (debug).
# - booster [default= gbtree ]Which booster to use. Can be gbtree, gblinear or dart; gbtree and dart use tree based models while gblinear uses linear fun
#--- AUX variables
# infolder <- "./outputs/"
# epochs            <- 150
# e                 <- 0.3
# gm                <- 0
# md                <- 6
# mcwt              <- 1
# ss                <- 1
# csbt              <- 1

Sys.setenv(CUDA="11.1") # Set global variable CUDA to 11.1 since pytorch and maybe other GPU-capable algorithms support only 11.1 version

# Set working directory to where the script is located
wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
#------------------ Requirements, dependencies and libraries
pkgs <- c('pacman','RMySQL','docopt','dplyr','xgboost', 'pacman', 'data.table','caret','ggplot2','tidyr','stringr')

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

"Mesothelial Cancer Cell Line Data Analysis (M.C.Li.D.A) Pipeline - Explants

Usage: explants.R --infolder=<folder> --epochs=<value> --max_depth=<value> --eta=<value> --gamma=<value> --colsample_bytree=<value> --min_child_weight=<value> --subsample=<value> --verbose=<value>

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

arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

#------------------ Functions
source("./pipeline/explants_fun.R")

#------------------ Load dataset and parameters into R environment 
infolder         <- arguments$infolder
epochs           <- as.numeric(arguments$epochs)
md               <- as.numeric(arguments$max_depth)
e                <- as.numeric(arguments$eta)           
g                <- as.numeric(arguments$gamma)
csbt             <- as.numeric(arguments$colsample_bytree)
mcwt             <- as.numeric(arguments$min_child_weight)
ss               <- as.numeric(arguments$subsample)   
verbose          <- as.numeric(arguments$verbose)

files         <- list(cfile = "./datasets/Explants/Clusters.csv",
                      sfile = "./datasets/Explants/Summary.xlsx" )



cluster_embed <- load_dataset(files$cfile)
sensitivity   <- load_dataset(files$sfile)
mname   <- "XGBoost_explants_logloss.model"
paths <- list(pdf       = "./outputs/pdf/", 
              models    = "./outputs/models/", 
              m.eval    = "./outputs/models/model_eval/", 
              r.objects = "./outputs/r-objects/"
)


#------------------ Connect to MySQL database and extract tables
details <- read.csv("../../mySQL/con", sep = " ")

c <- dbConnect(MySQL(), user=details$loginu, password=details$loginp, dbname=details$db, host=details$host) # Create connection to a local database server

dbListTables(c) # List tables
g <- dbGetQuery(c, "SELECT*FROM GainsEXP")
l <- dbGetQuery(c, "SELECT*FROM LossExp")
dc <- dbDisconnect(c)
#explants <- paste0(infolder, "MEDUSA50 DOMAINS - Original file.xlsx")

print(paste0("Project folder path set to: ", wd))
print(paste0("Extracting inputs from explants and saving sliced results into '", infolder,"'"))

gg  <- g[1:c(nrow(g)-4),] # 
gs <- g[c(nrow(g)-3):nrow(g),]

lg <- l[1:c(nrow(l)-4),]
ls <- l[c(nrow(l)-3):nrow(l),]

df <- merge(gg,lg)
df_deletions <- merge(gs,ls)
if(length(gg)+length(lg)-2 == length(df)){
  paste0("Succesefully merged the dataset into a new dataframe with a length of: ", length(df))
}

#------------------ Include evolutionary trajectory cluster embeds and IC50/logIC50 values
sensitivity           <- sensitivity$Summary
sensitivity           <- sensitivity[-9,]
colnames(sensitivity) <- sensitivity[1,]
sensitivity           <- sensitivity[-1,]
# inputs                <- sensitivity[,5:ncol(sensitivity)]
sensitivity           <- sensitivity[,c(1:4)]

explants_subsample <- sensitivity$Patient
cluster_embed <- cluster_embed[,-1]

cluster_embed <- cluster_embed[match(explants_subsample, cluster_embed[,1]),]
colnames(cluster_embed) <- c("sampleID", "Cluster")
cluster_embed$Cluster <- recode(cluster_embed$Cluster, "C1" = 1, "C2" = 2, "C3" = 3, "C4" = 4, "C5" = 5)

# explants_subsample %in% dataset
# cleaned <- dataset[match(explants_subsample, dataset[,1]),]

#------------------ Prepare labels
o <- c("inputs", "sample", "response")

# explants <- sliceDataset(inputs, o[1], genes)
# sampleID <- sliceDataset(sensitivity, o[2])

sampleID <- explants_subsample
labels <- list(response = as.numeric(sensitivity$response), 
               logIC50 = as.numeric(sensitivity$LogIC50), 
               IC50 = as.numeric(sensitivity$IC50))

#------------------- Prepare Inputs
match(sampleID, df$`medusa #`)
df <- df[match(sampleID, df$`medusa #`),]
inputs <- df[,-c(1,2)] # Slice uneccessary data from the original data frame


inputs <- toNumeric(inputs, o[1])
dataframe <- as.data.frame(cbind(sampleID, 
                                 labels,
                                 Evo_Clusters = cluster_embed[,-1], 
                                 inputs)) # Bind everything together to form the processed dataframe
pipe.data <- list(
  inputs    = inputs,
  labels    = labels,
  sid       = sampleID,
  dataframe = dataframe
)
saveRDS(pipe.data$dataframe,paths[1])
str(pipe.data$dataframe)
#------------------- Create a 70-30% Sample
set.seed(310)
qty <- 0.7
learning.set <- sampleDataset(qty,pipe.data$dataframe)
learning.set$test[,-c(1:4)]

#------------------ Create XGB matrices
train.matrix <- list(response = xgb.DMatrix(data = as.matrix(learning.set$train[,-c(1:4)]), 
                                            label = learning.set$train$response, missing = 99), 
                     logIC50  = xgb.DMatrix(data = as.matrix(learning.set$train[,-c(1:4)]), 
                                            label = learning.set$train$logIC50, missing = 99),
                     IC50     = xgb.DMatrix(data = as.matrix(learning.set$train[,-c(1:4)]), 
                                            label = learning.set$train$IC50, missing = 99)) # Convert into a list of
# dense matrices, set NA's as missing values to be imputed by xgboost

test.matrix  <- list(response = xgb.DMatrix(data = as.matrix(learning.set$test[,-c(1:4)]), 
                                            label = learning.set$test$response, missing = 99), 
                     logIC50  = xgb.DMatrix(data = as.matrix(learning.set$test[,-c(1:4)]), 
                                            label = learning.set$test$logIC50, missing = 99),
                     IC50     = xgb.DMatrix(data = as.matrix(learning.set$test[,-c(1:4)]), 
                                            label = learning.set$test$logIC50, missing = 99))

#------------------ Set the parameters
params    <- list(response = list(booster          = "gbtree",
                                  eval_metric      = "mlogloss",
                                  objective        = "multi:softmax",
                                  num_class        = 5,
                                  eta              = e,
                                  gamma            = gm,
                                  min_child_weight = mcwt,
                                  max_depth        = md,
                                  colsample_bytree = csbt,
                                  subsample        = ss),
                  logIC50  = list(booster          = "gbtree",
                                  eta              = e),
                  IC50     = list(booster          = "gbtree",
                                  eta              = e))

save_name <- mname
seed      <- 310
watchlist <- list(response = list(train = train.matrix$response, 
                                  test = test.matrix$response),
                  logIC50  = list(train = train.matrix$logIC50, 
                                  test = test.matrix$logIC50),
                  IC50     = list(train = train.matrix$IC50, 
                                  test = test.matrix$IC50))

#------------------ Run Extreme Gradient Boosting algorithm (XGBOOST)
medusa <- list(response.model = createModel(data = train.matrix$response, 
                                            params = params$response, 
                                            watchlist = watchlist$response, 
                                            epochs = epochs, 
                                            save_name = save_name, 
                                            seed = seed),
               logIC50.model  = createModel(data = train.matrix$logIC50, 
                                            params = params$logIC50, 
                                            watchlist = watchlist$logIC50, 
                                            epochs = epochs, 
                                            save_name = save_name, 
                                            seed = seed),
               IC50           = createModel(data = train.matrix$IC50, 
                                            params = params$IC50, 
                                            watchlist = watchlist$IC50, 
                                            epochs = epochs, 
                                            save_name = save_name, 
                                            seed = seed))

e.log.list <- list(r.log = data.frame(medusa$response.model$evaluation_log),
                   lic50.log = data.frame(medusa$logIC50.model$evaluation_log),
                   ic50.log = data.frame(medusa$IC50$evaluation_log))

# Extract the xgboost training evaluation log
imatrices <- list(r.imp     = xgb.importance(colnames(learning.set$train[,-c(1:4)]), medusa$response.model), # Extract the matrix of feature importance
                  lic50.imp = xgb.importance(colnames(learning.set$train[,-c(1:4)]), medusa$logIC50.model),
                  ic50.imp  = xgb.importance(colnames(learning.set$train[,-c(1:4)]), medusa$IC50))

xgb.plot.shap(data = as.matrix(learning.set$train[,-c(1:4)]), model = medusa$response.model, top_n = 10, plot_NA = TRUE)
xgb.plot.shap(data = as.matrix(learning.set$train[,-c(1:4)]), model = medusa$logIC50.model, top_n = 10, plot_NA = FALSE)
xgb.plot.shap(data = as.matrix(learning.set$train[,-c(1:4)]), model = medusa$IC50, top_n = 10, plot_NA = FALSE)

#------------------ Evaluation plot
pdf(paste0(paths[1],"summary.pdf"), height = 11.69, width = 8.27)
r.plots     <- p(e.log.list$r.log, medusa$response.model, paths, as.matrix(learning.set$train[,-c(1:4)]), imatrices$r.imp)
lic50.plots <- p(e.log.list$lic50.log, medusa$logIC50.model, paths, as.matrix(learning.set$train[,-c(1:4)]), imatrices$lic50.imp)
ic5.plots <- p(e.log.list$ic50, medusa$IC50, paths, as.matrix(learning.set$train[,-c(1:4)]), imatrices$ic50.imp)
dev.off()
#------------------ Model Predictions
print(paste0("Predicting on explants via XGBoost with inputs from ", infolder))
medusa.predict <- list(response = predict(medusa$response.model, newdata = as.matrix(learning.set$test[,-c(1:4)]), missing = 99),
                       LogIC50  = predict(medusa$logIC50.model, newdata = as.matrix(learning.set$test[,-c(1:4)]), missing = 99),
                       IC50     = predict(medusa$IC50, newdata = as.matrix(learning.set$test[,-c(1:4)]), missing = 99))
medusa.predict$response
learning.set$test$response
medusa.predict$LogIC50
learning.set$test$logIC50

medusa.predict$IC50
learning.set$test$IC50
metrics <- calcMetrics(medusa.predict$response, learning.set$test$response)

out <- data.frame(cbind(sampleID           = learning.set$test$sampleID, 
                        response           = learning.set$test$response,
                        predicted.response = medusa.predict$response,
                        LogIC50            = learning.set$test$logIC50,
                        predicted.LogIC50  = round(medusa.predict$LogIC50,2),
                        IC50               = learning.set$test$IC50,
                        predicted.IC50     = round(medusa.predict$IC50,2)
))

write.csv(out, paste0(infolder,"output.csv"))

#------------------ Save the model and other processed stuff
model.names <- c("medusa.response", "medusa.LogIC50", "medusa.IC50")
i <- 1
for(models in medusa){
  xgb.save(models, paste0(paths$models, model.names[i]))
  i <- i + 1
}

print(paste0("Files were succesefully saved as '", model.names[1],"', '",model.names[2], "' and '", model.names[3],"'." ))


