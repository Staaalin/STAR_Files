#!/bin/csh


set Energy = 62.0
# set midname = "/star/data01/pwg/svianping/output/output_"
set midname = "/star/data01/pwg/svianping/HADD/HADD_T_"
set outname = "~/Result/K0SXiMix/HADDr_"
set StartFileIndex = 0
set EndFileIndex = 341
set OutputFileIndex = 1
set A_PDG = 310
set B_PDG = 3312

root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,1\)
