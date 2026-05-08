#!/bin/csh


set midname = "/star/data01/pwg/svianping/output/output_"
# set midname = "/star/data01/pwg/svianping/HADD/HADD_T_"
# set midname = "~/Result/Cor_T_"

set outname = "/star/data01/pwg/svianping/HADD/HADD_"
# set outname = "~/Result/Cor_"
set StartFileIndex = 151
set EndFileIndex = 250
set OutputFileIndex = 52
set A_PDG = \-321
set B_PDG = 3122

# root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,0\)
# root -b MM.C\(\"$midname\",\"dAu_200_21\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,0\)
# root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,1\)
# root4star -b CheckParent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\"\)
# root4star -b MixEventTest.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,0\)

./Side \
"FileList" \
"dAu_200_21" \
151 250 52 \
"/star/data01/pwg/svianping/HADD/HADD_" \
211 2212 \-211 0 0
