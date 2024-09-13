#!/bin/csh


# set midname = "/star/data01/pwg/svianping/output/output_"
set midname = "/star/data01/pwg/svianping/HADD/HADD_T_"
# set midname = "~/Result/Cor_T_"

set outname = "/star/data01/pwg/svianping/HADD/HADD_"
# set outname = "~/Result/Cor_"
set StartFileIndex = 0
set EndFileIndex = 51
set OutputFileIndex = 52
set A_PDG = 321
set B_PDG = 3122

# root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,0\)
root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,1\)
