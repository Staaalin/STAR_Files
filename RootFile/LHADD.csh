#!/bin/csh

cd /star/data01/pwg/svianping/output/

set numFiles = 'ls -l output_1000*.root | wc -l'

echo $numFiles

# set i = 0
# set j = 0
# while ($i <= $numFiles)

#     source /star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh
#     set ARM = " > HADD_"
#     root4star -q -b HADDr_xml.C\(\"$ObvInputName\",\"$ObvOutputName\",$i,$FilesPerJob,$FileStart,$FileEnd\)$ARM$i.log

#     @ i = $i + 1
# end
