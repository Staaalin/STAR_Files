#!/bin/csh


set Energy = 62.0
# set midname = "/star/data01/pwg/svianping/output/output_"
set midname = "~/Result/KLaMix/HADDr_"
set outname = "~/Result/KLaMix/HADDr_"
set StartFileIndex = 0
set EndFileIndex = 0
set OutputFileIndex = 1
set A_PDG = 321
set B_PDG = 3122

root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,1\)
