#--- AUX
infolder <- "./outputs/r-objects/"

pkgs <- c('pacman','dplyr','docopt','data.table', 'caTools')

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

#------------------Functions
source("./pipeline/patients_fun.R")

#------------------ File paths to maintain the file system integrity
mname   <- "MiST-XGBoost.model"
paths <- list(pdf = "./outputs/pdf/")
dfile  <- paste0(infolder, "dataframe.RDS")

#------------------ Read the saved R object
dataframe <- readRDS(dfile)

#------------------ Sample the dataset (70%-30%)
set.seed(312)
qty <- 0.7
learning.set <- sampleDataset(qty,dataframe)
learning.set$test

#------------------ Run LogitBoost algorithm
?LogitBoost
t <- LogitBoost(dataframe[,-1],  learning.set$train$labels, nIter = 100)
t
dataframe[,-1]
### EXAMPLE ###

data(iris)
data = iris[,-5]
label = iris[, 5]
data
label
model <- LogitBoost(data, label, nIter = 20)
model
