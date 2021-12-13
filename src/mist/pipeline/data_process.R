#-------------------------------------------------------------------------------
# Evolutionary Trajectory-Drug Sensitivity via ENSEMBLE Learning Pipeline
#-------------------------------------------------------------------------------
# 0   lack		
# -1	loss -> WE expect this to be  Associated with sensitivity. 	Bridge analogy for  Synthetic lethality as a response. 
# -2	deletion		
# 1	  gain		
# 2	  amplifications		
# 3	  upd	-> Unparental  Disomy  (mimics deletion/amplification)	==> SET THIS TO NA and we reduce the number of dimensions to 4
# 3 - Complete Response (CR): Disappearance of all target lesions
# 2 - Partial Response (PR): At least a 30% decrease in the sum of the LD of target lesions, taking as reference the baseline sum LD
# 0 - Progressive Disease (PD): At least a 20% increase in the sum of the LD of target lesions, taking as reference the smallest sum LD recorded since the treatment started or the appearance of one or more new lesions
# 1 - Stable Disease (SD): Neither sufficient shrinkage to qualify for PR nor sufficient increase to qualify for PD, taking as reference the smallest sum LD since the treatment started

#--- AUX
# infolder <- "./outputs/"
# path     <- "/home/jr429/Documents/UoL/Bioinformatics_MSc/IndependentResearchProject/ET-DS-via-ENSEMBLE-IRP/src/mist-meso-pipe/"
# setwd(path)

Sys.setenv(CUDA="11.1")

#------------------ Requirements, dependencies and libraries
pkgs <- c('pacman','docopt','dplyr','stringr', 'data.table', 'RMySQL','readxl')


suppressMessages(if (!require("BiocManager", character.only = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
} else {
  ipkgs <- sapply(pkgs, function(...) require(..., character.only = TRUE))
  if (any(!ipkgs)) {
    BiocManager::install(pkgs[!ipkgs])
    install.packages(pkgs[!ipkgs])
  } else {
    message("\n\nCool! Your machine has everything is needed.\n\n")
  }
})


print("Loading required packages...")
library(pacman)
pacman::p_load(pkgs, install = TRUE, character.only = TRUE)
pacman::p_loaded()

"Mesothelial Cancer Cell Line Data Analysis (M.C.Li.D.A) Pipeline - Data preparation

Usage: data_process.R --infolder=<folder> --path=<folder>

Options:
  -h --help                  Show this screen
  --infolder=<folder>        Folder where with the are placed
  --path=<folder>            Path to the pipeline project folder (i.e. if you downloaded into downloads its '/home/user/Downloads/mist-meso-pipe/'
  --version                  Pipeline Pre-release version
"-> doc

arguments <- docopt(doc, quoted_args = TRUE, help = TRUE)
print(arguments)

#------------------ Functions
source("./pipeline/funtions.R")

#------------------ Load dataset and parameters into R environment 
infolder <- arguments$infolder
path     <- arguments$path
wd       <- setwd(path)
names    <- list(iname = "inputs.RDS", 
                 lname = "labels.RDS",
                 sname = "sid.RDS",
                 cname = "dataframe.RDS",
                 mname = "metadata.RDS")


print(paste0("Project folder path set to: ", wd))
print(paste0("Preparing patients dataset with ", infolder, " as the folder where slices are generated!"))

pfile <- read_excel("./datasets/Patients/DATAFRAMES MIST1.xlsx")
CNV   <- pfile[,-c(79:ncol(pfile))]

#------------------ Label the training data (and for convenience, slice data set into inputs, sample ID and labels)
mv        <- apply(CNV[,-c(1:5)], 2, calcMissv)        # Calculate % missing data
geneset   <- colnames(CNV[,grep("^[A-Z]([A-Z]|\\d)+", colnames(CNV))]) # Extract genes from the corresponding column names
paste0("Summary of missing data: ", round(mean(mv),2), " Standard deviation: ", round(sd(mv),2))

o <- c("inputs", "sample", "response")
inputs    <- sliceDataset(CNV, o[1]) 
sample    <- sliceDataset(CNV, o[2])
labels    <- sliceDataset(CNV, o[3])

inputs    <- toNumeric(as.data.frame(inputs), o[1])
labels    <- toNumeric(labels, o[3])

data      <- list(inputs    = inputs,  # List of cleaned data
                  labels    = labels, 
                  sample    = sample,  
                  dataframe = cbind(sample,labels,inputs), # This is the final processed data frame to be used
                  meta.data = list(gene.set = geneset, missing.values = mv)) # Some extra informations

#------------------ Some auxiliary processing
#inputs$type <- sample # Add a seperate column for test type
#inputs$type <- str_extract(inputs$type, "\\w$") # Extract test type as the last character in sample name into "type" column
#sample      <- str_extract(sample, "\\w+\\d+") # Extract sample names into "sample" vector
#inputs$type <- type[inputs$type] # Convert test types into numerical values where Test is represented by 0 and Placebo as 1

#------------------ Transform factor into dummy variable (NO NEED)
#install.packages("fastDummies")
#library("fastDummies")
#x <- dummy_cols(x, remove_first_dummy = TRUE)
#str(x)

#------------------ Save the processed table as an R object to outputs folder
files <- list(ifile=paste0(infolder, "r-objects/", names$iname),
              lfile=paste0(infolder, "r-objects/", names$lname),
              sfile=paste0(infolder, "r-objects/", names$sname),
              cfile=paste0(infolder, "r-objects/", names$cname),
              mfile=paste0(infolder, "r-objects/", names$mname))
i <- 1
for(n in files){
  print(i)
  saveRDS(data[[i]], files[[i]])
  i <- i + 1
}

if(file.exists(files[[1]]) && file.exists(files[[2]]) && file.exists(files[[3]])){
  print(paste0("Files were succesefully saved as '", files[1],"', '", files[2],"', ",files[3], "' and ","'",files[4],"'"))
}
