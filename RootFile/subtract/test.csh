#!/bin/csh


set Energy = 62.0
set midname = "/star/data01/pwg/svianping/output/output_"
set StartFileIndex = 500
set EndFileIndex = 500
set OutputFileIndex = 0
set A_PDG = 321
set B_PDG = 3312

# root4star -b Subtract.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,$Energy\)
root4star -b MixEvent\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,$A_PDG,$B_PDG\)