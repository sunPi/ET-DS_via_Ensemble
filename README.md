# ET-DS_via_EnsembleET-DS_via_Ensemble
This repository is to track versioning of the research pipeline for exploring evolutionary trajectory-drug sensitivity interactions through ensemble learning method. Currently, XGBoost algorithm is implemented. The algorithm is used to trai and evaluate the model. Following this, it is used to predict. In parallel with training, knowledge that the machine learns is extracted and explained via SHAP values (Explained here: https://github.com/slundberg/shap). The scores measure feature importance across all classs (global explainability) and for each individual class of each feacture (local explainability).

The pipeline:

There are three scripts, connected together via a shell script. The scripts either extract and process data, train the model and extract knowledge and ultimatively predict the responses. The pipeline is created with a preserved file structure and it should be cloned as a whole to work properly.

1) Data preprocessing
The script extracts relevant gene features from a MySQL database (currently implemented only for local db). It then slices the relevant information, cleans any uneccessary data and encodes categorical variables appropriately. The output of the script are several tables mainly inputs, labels, patienID, metadata and a data frame consisting of all previously mentioned tables. They are saved as .RDS files in the output/r-objects folder.

2) Model Training/Evaluation
Training:

XGBoost's R implementation is used to train and evaluate the model (via k-fold cross-validation). This script loads saved RDS objects from script 1 and creates a learning set list, splitting the data into two sets: .7 goes to a training set, .3 goes to the testing set, all row wise. Following that, sets are converted into dense DMatrices, that XGBoost supports. The value "99" is set as the missing values. Watchlist is created for both sets, that tracks model performance during training and evaluation. A list of training hyper parameters is created, starting of with the default values (i.e. eta = 0.3 etc.), and nrounds is set to 150. Then, the model is trained using the dense training matrix, parameters list and watchlist.

Evaluation:

Since this is a tree-based algorithm, evaluation is necessary to reduce the chance of overfitting. A native k-fold cross-validation is implemented in the XGBoost package, that is used after the first trained model. Currently, the evaluation function trains the model until it has not improved yet in the next 20 rounds. It then plots the evaluation log, saves it to a pdf file. It then finds the training iteration, where test-error was the lowest and updates the hyperparameter "nrounds". The model is retrained on new parameters.

In the future, there is plan to repeat the process for other hyperparameters as well.

The model is saved as an "XGBoost model" file, which is a package-native format. It can be loaded only via the same pagackage load function. Other outputs include graphs of basic/trained model evaluation logs, model complexity measures, SHAP summary and local/global feature importance, that are saved to the outputs/pdf folder.

3) Model Testing
Output of the previous script, the model is then loaded via the package functions. A custom function to calculate metrics such as accuracy, precision, recall and F1 score was written to make the pipeline more compact. The metrics, together with the predictions are joined into a data table, that is then saved as a comma seperated value file.

How to use the pipeline:

First, clone the **whole** github repository into a local folder then open a linux terminal and run the run_ETDS.sh script.
