#!/bin/csh


set Energy = 62.0
# set midname = "/star/data01/pwg/svianping/output/output_"
set midname = "/star/data01/pwg/svianping/HADD/HADD_"
set outname = "/star/data01/pwg/svianping/HADD/HADDr_"
set StartFileIndex = 0
set EndFileIndex = 27
set OutputFileIndex = 0
set A_PDG = 321
set B_PDG = 3312

<<<<<<< HEAD
root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG,0\)
=======
root4star -b MixEvent.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,\"$outname\",$A_PDG,$B_PDG\)
>>>>>>> parent of 13ac31b2 (upcode)
