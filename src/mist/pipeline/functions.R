#------------------ Data Preparation Functions
load_dataset <- function(dataset) {
  print("-------------------------------")
  print("Loading the dataset into R...")
  if(grepl("\\.RDS$", dataset)){
    dataset   <- readRDS(dataset)
    return(as.data.frame(dataset))
  } else if(grepl("\\.xlsx$", dataset)){
    print("Please implement the read excel function!")
    #dataset   <- fread(dataset)
    #return(as.data.frame(dataset))
  } else{
    dataset   <- readRDS(dataset)
    return(as.data.frame(dataset))
  }
}
extractSQL <- function(statements, details){
  c <- dbConnect(MySQL(), user=details$loginu, password=details$loginp, dbname=details$db, host=details$host)
  dbListTables(c) # List tables
  table <- dbGetQuery(c, statements)
  dc <- dbDisconnect(c)
  return(table)
}
calcMissv <- function(x){
  sum(is.na(x))/length(x)*100
}
sliceDataset <- function(dataset, type){
  if(type == "inputs"){
    toInputs(dataset, genes)
  } else if(type == "sample"){
    toSample(dataset)
  }
  else if(type == "response"){
    toResponse(dataset)
  }
}
toInputs <- function(dataset=dataset, genes=genes, unique=unique){
  #colnames(dataset) <- colnames(dataset) %>% make.unique() # Make duplicated column names into unique names
  inputs <- dataset[,grep("^[A-Z]([A-Z]|\\d)+", colnames(CNV))]
  return(inputs) 
}
toSample <- function(dataset=dataset) {
  sample <- dataset$`Patient ID`
  return(sample)
}
toResponse <- function(dataset=dataset){
  response <- dataset$`response category`
}
toNumeric <- function(dataset,option){
  #dataset[is.na(dataset)] <- NULL
  #dataset[dataset == ""] <- NA # Insert NA where cells are empty
  #sum(is.na(dataset))/length(dataset)*100
  if(option == "inputs"){
    dataset <- as.data.frame(dataset)
    for (i in colnames(dataset)) {
      dataset[,i] <- dplyr::recode(dataset[,i], "deletion" = -2, "loss" = -1, "wild type" = 0, "gain" = 1, "amplification" = 2, .default = 99, .missing = 99)
    }
  }
  else if(option == "binary"){
    for (i in names(dataset)) {
      dataset[,i] <- dplyr::recode(dataset[,i], "deletion" = 0, "loss" = 1, "wild type" = 0, "gain" = 0, "amplification" = 0, .default = 99, .missing = 99)
    }
  }
  else if(option == "response"){
    dataset <- dplyr::recode(dataset, "SD" = 1, "PD" = 0, "PR" = 1, "CR" = 1, .default = NULL, .missing = NULL)
  }
  return(dataset)
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
createModel <- function(data, inputs, labels, params, watchlist, epochs, save_name, verbose, seed){
  if(!missing(data)){
    xgb.train(data      = data, 
              params    = params, 
              watchlist = watchlist,
              nrounds   = epochs,
              save_name = save_name,
              verbose   = verbose,
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
    confusion.matrix = confusionMatrix(as.factor(predictions), as.factor(labels)),
    accuracy         = round((sum((predictions == labels), na.rm = TRUE))/length(labels), 2),
    precision        = precision(as.factor(predictions), as.factor(labels)),
    recall           = recall(as.factor(predictions), as.factor(labels)),
    F1               = F_meas(as.factor(predictions), as.factor(labels))
  )
  return(metrics)
}
p <- function(e, model, paths, train.matrix, imatrix){
  plots <- list(logloss = list(plot(e$iter, e$train_mlogloss, col = 'blue', xlab = "No. of iterations", ylab = "Log-loss of the likelihood function"), 
                               lines(e$iter, e$test_mlogloss, col = 'red'),
                               title("Training (blue) + Validation (red) logloss change")),
                shap    = xgb.plot.shap(data = train.matrix, model = model, top_n = 10, plot_NA = TRUE), # Extract top 10 features based on their SHAP importance values
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