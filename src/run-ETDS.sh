#!/bin/bash

# ---------------------
# Evolutionary Trajectory-Drug Sensitivity via ENSEMBLE Learning Pipeline for Mesothelioma
# ---------------------

#--------- Run Mist Pipeline
#sh /home/jr429/Documents/UoL/Bioinformatics_MSc/IndependentResearchProject/ET-DS_via_Ensemble/src/mist/run-mist.sh

#--------- Run Explants Pipeline

#cd explants
#sh run-explants.sh

parallel  ::: 'cd mist && ./run-mist.sh' 'cd explants && ./run-explants.sh'
