#!/bin/bash
# ---------------------
# Evolutionary Trajectory-Drug Sensitivity via ENSEMBLE Learning Pipeline for Mesothelioma
# ---------------------

# ---------------------
# Format Patient MySQL Tables into MiST Pipeline XGBoost-ready Format
# ---------------------
echo Running explants pipeline...

infolder="./outputs/"
epochs=150
max_depth=6
eta=0.3
gamma=0
colsample_bytree=1
min_child_weight=1
subsample=1
verbose=0

path=$(dirname "$(realpath IRP-pipeline.sh)")
echo $path

Rscript ./pipeline/explants.R --path="${path}" --infolder="${infolder}" --epochs="${epochs}" --max_depth=${max_depth} --eta=${eta} --gamma=${gamma} --colsample_bytree=${colsample_bytree} --min_child_weight=${min_child_weight} --subsample=${subsample} --verbose=${verbose}
