#------------------ Global
requirements <- function(pkgs){
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
  return(pacman::p_loaded())
}
#------------------ Data Preparation Explants
load_dataset <- function(dataset) {
  print("-------------------------------")
  print("Loading the dataset into R...")
  if(!grepl("\\.xlsx$",dataset)){
    dataset   <- fread(dataset)
    return(as.data.frame(dataset))
  }
  
  if(grepl("\\.xlsx$", dataset)){
    pacman::p_load(readxl)
    es <- excel_sheets(dataset)
    if(es > 1){
      i <- 1
      for(s in es){
        dataset[i] <- list(as.data.frame(read_excel(dataset, es[i])))
        i <- i+1
      }
      names(dataset) <- es
    }  else{
      dataset <- read_excel(dataset)
    }
    pacman::p_unload(readxl)
    return(dataset)
  }
}

sliceDataset <- function(dataset, type, genes){
  if(type == "inputs"){
    toInputs(dataset, genes)
  } else if(type == "sample"){
    toSample(dataset)
  }
  else if(type == "response"){
    toResponse(dataset)
  }
}
toInputs <- function(dataset=dataset, genes=genes){
  colnames(dataset) <- colnames(dataset) %>% make.unique() # Make duplicated column names into unique names
  inputs <- dataset[-1,]
  inputs <- inputs[order(inputs[,2]),]
  inputs <- inputs[,-c(1:2)]
  #inputs <- dataset %>% select(diff) # Make inputs based on a list of genes
  return(inputs) 
}
toSample <- function(dataset=dataset) {
  sample <- dataset[-1,2]
  return(sample)
}
toResponse <- function(dataset=dataset){
  response <- dataset[,4]
}
toNumeric <- function(dataset,option){
  #dataset[dataset == ""] <- NA # Insert NA where cells are empty
  if(option == "inputs"){
    for (i in names(dataset)) {
      dataset[,i] <- recode(dataset[,i], "Sd" = -2, "Cd" = -1, "wild-type" = 0, "C" = 1, "S" = 2, .default = 0, .missing = NULL)
    }
  }
  else if(option == "response"){
    dataset <- recode(dataset, "PD" = 0, "SD" = 1, "PR" = 1, "CR" = 1, .default = NULL, .missing = NULL)
  }
  return(dataset)
}
cutTail <- function(lists){
  i <- 1
  for(list in lists){
    #cut <- as.numeric(rownames(tail(lists[[i]], 4)))
    df[[i]] <- as.data.frame(lists[[i]][-c(48:51),])
    i <- i+1
  }
  return(df)
}
sampleDataset <- function(percent, y){
  if(is.data.frame(y)){
    sqty      <- round(percent*nrow(y))
    indx      <- sample(nrow(y),sqty)
    sample    <- y[indx,]
    remainder <- y[-indx,]
    return(list(train = sample, test = remainder))
  }
  else if(is.vector(y)){
    sqty   <- round(percent*length(y))
    sample <- y[sample(nrow(y),sqty),]
    return(sample)
  }
}

#------------------ Machine Learning Functions
createModel <- function(data, inputs, labels, params, watchlist, epochs, save_name, seed){
  if(!missing(data)){
    xgb.train(data      = data, 
              params    = params, 
              watchlist = watchlist,
              nrounds   = epochs,
              save_name = save_name,
              set.seed(seed))
  }
  else if(!missing(inputs) && !missing(labels)){
    xgboost(data = inputs,
            label = labels,
            params = params, 
            watchlist,
            nrounds = epochs,
            save_name = save_name,
            set.seed(311))
  }
}
calcMetrics <- function(predictions, labels){
  metrics = list(
    accuracy         = round((sum((predictions == labels), na.rm = TRUE))/length(labels), 2),
    precision        = precision(as.factor(predictions), as.factor(labels)),
    recall           = recall(as.factor(predictions), as.factor(labels)),
    F1               = F_meas(as.factor(predictions), as.factor(labels)),
    if(0 %in% labels && 1 %in% labels){
      confusion.matrix = confusionMatrix(as.factor(predictions), as.factor(labels))
    }
  )
  return(metrics)
}
p <- function(e, model, paths, train.matrix, imatrix){
  plots <- list(logloss = list(plot(e$iter, e$train_mlogloss, col = 'blue', xlab = "No. of iterations", ylab = "Log-loss of the likelihood function"), 
                               lines(e$iter, e$test_mlogloss, col = 'red'),
                               title("Training (blue) + Validation (red) logloss change")),
                shap    = xgb.plot.shap(data = train.matrix, model = model, top_n = 10), # Extract top 10 features based on their SHAP importance values
                complex = xgb.ggplot.deepness(model),                                                              # Complexity measurement
                ft.cont = xgb.ggplot.importance(imatrix),                                                          # Feature contribution measurement
                sh.summ = xgb.ggplot.shap.summary(train.matrix, model = model, top_n = 20))
  print(plots)
}
pCV <- function(paths, e.cv, e.cv.re){
  if(!missing(e.cv)){
    pdf(paste0(paths, "eval.pdf"), height = 11.69, width = 8.27)
    scatter.smooth(e.cv, xlab = "No. of Rounds", 
                   ylab = "Change in mLogLoss - CV") +
      lines(e.cv, col = "red")
    dev.off()
  }
  else if(!missing(e.cv.re)){
    pdf(paste0(paths, "eval_re.pdf"), height = 11.69, width = 8.27)
    scatter.smooth(e.cv.re, xlab = "No. of Rounds", ylab = "Change in mLogLoss")+
      lines(e.cv.re$iter, e.cv.re$test_mlogloss, col = 'red')
    dev.off()
  }
}

crossVal <- function(model.cv){
  lv.iter  <- min(model.cv$evaluation_log$test_mlogloss_mean)
  lv.index <- which(model.cv$evaluation_log$test_mlogloss_mean == lv.iter)
  value    <- model.cv$evaluation_log$iter[lv.index] # This variable indicates the lowest point of error lost in the testing mean
  return(value)
}