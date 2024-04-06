#!/bin/csh

echo "Please enter which location:"
echo "SCHEME 1: /star/data01/pwg/svianping/output/output_*.root"
echo "SCHEME 2: /star/data01/pwg/svianping/HADD/HADD_*.root"
echo "SCHEME 3: /star/data01/pwg/svianping/HADD/HADDrA_*.root"
set InputNameIndex = "$<"

echo "Please enter merge how much .root into ONE:"
set FilesPerJob = "$<"

echo "Please enter from which file:"
set FileStart = "$<"

echo "Please enter to which file:"
set FileEnd = "$<"

if ($InputNameIndex == 1) then
    set ObvInputName = "/star/data01/pwg/svianping/output/output_"
    set ObvOutputName = "/star/data01/pwg/svianping/HADD/HADD_"
    set OutputURL = "/star/data01/pwg/svianping/HADD/"
else if ($InputNameIndex == 2) then
    set ObvInputName = "/star/data01/pwg/svianping/HADD/HADD_"
    set ObvOutputName = "/star/data01/pwg/svianping/HADD/HADDrA_"
    set OutputURL = "/star/data01/pwg/svianping/HADD/"
else if ($InputNameIndex == 3) then
    set ObvInputName = "/star/data01/pwg/svianping/HADD/HADDrA_"
    set ObvOutputName = "/star/u/svianping/STAR_Files/RootFile/HADDrB_"
    set OutputURL = "/star/u/svianping/STAR_Files/RootFile/"
else
    echo "Error INVALID location!"
    exit
endif

set MainDir=`pwd`


set numFiles = 1
@ numFiles = ( $FileEnd - $FileStart ) / $FilesPerJob

set i = 0
set j = 0
while ($i <= $numFiles)

    echo "source /star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh"
    set ARM = " > HADD_"
    root4star -q -b HADDr_xml.C\(\"$InputName\","$OutputName",$i,$FilesPerJob,$FileStart,$FileEnd\)$ARM$i.log

    @ i = $i + 1
end
