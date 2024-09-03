#!/bin/csh


set midname = "/star/data01/pwg/svianping/HADD/HADD_T_"
set outname = "~/Result/Cor_"
set StartFileIndex = 0
set EndFileIndex = 341
set OutputFileIndex = 1
set A_PDG = 321
set B_PDG = 3122

root4star -b MixEvent_Test\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,0\)
