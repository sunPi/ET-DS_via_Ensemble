#--- AUX
infolder <- "./outputs/r-objects/"

wd     <- setwd("~/Documents/UoL/Bioinformatics MSc/IndependentResearchProject/ET-DS-via-ENSEMBLE-IRP/")

pkgs <- c('pacman','dplyr','docopt','data.table', 'fastAdaboost', 'caret')

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

#------------------ Run fast Adaboost algorithm
learning.set$train$labels <- as.factor(learning.set$train$labels)
real_adaboost(tr[2]~., data , 10,)

learning.set$train[,3]

fakedata <- data.frame( X=c(rnorm(100,0,1),rnorm(100,1,1)), Y=c(rep(0,100),rep(1,100) ) )
fakedata$Y <- factor(fakedata$Y)
test_adaboost <- real_adaboost(Y~X, data=fakedata,10)
