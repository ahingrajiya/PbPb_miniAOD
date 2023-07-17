#!/bin/bash

echo "Setup CMSSW (Import ROOT version)"
cd /afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/
eval `scramv1 runtime -sh`
cd /afs/cern.ch/user/a/ahingraj/private/analysis/CMSSW_12_5_0/src/myProcesses/jettrackcorrelation_analyzer/PbPb_code
echo "Submit skim jobs at "
echo PWD: $PWD

root -l -b -q "PbPb_studies.cpp(\"$1\",\"$2\")"