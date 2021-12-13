#!/bin/bash
# ---------------------
# Evolutionary Trajectory-Drug Sensitivity via ENSEMBLE Learning Pipeline for Mesothelioma
# ---------------------

# Load Global Arguments
#jobname="Job-testing-the-pipeline"
#infolder="./datasets/"
#outfolder="./outputs/"

# ---------------------
# Extract Tables From MySQL Database
# ---------------------
#echo Running MySQL Extract Tables script...
#Rscript ./pipeline/MySQL-to-R.R --jobname=${jobname} --infolder=${infolder}  --outfolder="${outfolder}"

# ---------------------
# Format Patient MySQL Tables into MiST Pipeline XGBoost-ready Format
# ---------------------
echo Running data preparation for patients...

infolder="./outputs/"
path=$(dirname "$(realpath IRP-pipeline.sh)")
echo $path

Rscript ./pipeline/data_process.R --infolder="${infolder}" --path="${path}"

# ---------------------
# Reduce dimensions using several approaches
# ---------------------
# echo Reducing dimensions...
#
# infolder="./outputs/r-objects/"
#
# Rscript ./pipeline/MiST_PCA_patients.R --infolder="${infolder}"

# ---------------------
# Run XGBoost train on Patients
# ---------------------
echo Training the model...

infolder="./outputs/r-objects/"
epochs=150
max_depth=6
eta=0.3
gamma=0
colsample_bytree=1
min_child_weight=1
subsample=1
verbose=0

Rscript ./pipeline/train.R  --infolder="${infolder}" --epochs="${epochs}" --max_depth=${max_depth} --eta=${eta} --gamma=${gamma} --colsample_bytree=${colsample_bytree} --min_child_weight=${min_child_weight} --subsample=${subsample} --verbose=${verbose}

# ---------------------
# Run XGBoost predict on Patients
# ---------------------
infolder="./outputs/"
echo Predicting patient responses...

Rscript ./pipeline/test.R --infolder="${infolder}"

# ---------------------
# Format Explant MySQL Tables into XGBoost-ready format
# ---------------------

#Rscript ./pipeline/MiST-data-preparation_explants.R --infolder="${infolder}"


# ---------------------
# Run XGBoost predict on Explants
# ---------------------

#Rscript ./pipeline/CLI_MiST-XGBoost_predict.R --infolder="${infolder}"


# ---------------------
# Running SGD, XGBoost and Boruta in parallel
# ---------------------
#inputs="./outputs/minst-inputs.RDS"
#labels="./outputs/minst-labels.RDS"
#epochs=180
#dataset2="./datasets/ds.csv"
#column_targets="YieldA, YieldP"


#Rscript ./R-scripts/CLI_RandomDS-SGD.R --dataset="${dataset2}" --infolder="${infolder}" --column_targets="${column_targets}"
