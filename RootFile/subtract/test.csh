#!/bin/csh

set LambdaFitMass = 1.1161
set LambdaFitMassSigma = 0.002
set LambdaBarFitMass = 1.1161
set LambdaBarFitMassSigma = 0.002
set XiFitMass = 1.3223
set XiFitMassSigma = 0.0024
set XiBarFitMass = 1.3223
set XiBarFitMassSigma = 0.0024
set OmegaFitMass = 1.6725
set OmegaFitMassSigma = 0.0029
set OmegaBarFitMass = 1.6727
set OmegaBarFitMassSigma = 0.0024


set midname = "/star/data01/pwg/svianping/output/output_"
set StartFileIndex = 500
set EndFileIndex = 500
set OutputFileIndex = 0

root4star -b Subtract.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,$LambdaFitMass,$LambdaFitMassSigma,$XiFitMass,$XiFitMassSigma,$OmegaFitMass,$OmegaFitMassSigma,$LambdaBarFitMass,$LambdaBarFitMassSigma,$XiBarFitMass,$XiBarFitMassSigma,$OmegaBarFitMass,$OmegaBarFitMassSigma\)