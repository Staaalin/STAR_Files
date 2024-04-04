#!/bin/csh

# set InputName = 
# set FilesPerJob = 400
# set FileStart = 1
# set FileEnd = 66389

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
    set InputName = "output_"
    set OutputName = "HADD_"
else if ($InputNameIndex == 2) then
    set ObvInputName = "/star/data01/pwg/svianping/HADD/HADD_"
    set ObvOutputName = "/star/data01/pwg/svianping/HADD/HADDrA_"
    set InputName = "HADD_"
    set OutputName = "HADDrA_"
else if ($InputNameIndex == 3) then
    set ObvInputName = "/star/data01/pwg/svianping/HADD/HADDrA_"
    set ObvOutputName = "/star/u/svianping/STAR_Files/RootFile/HADDrB_"
    set InputName = "HADDrA_"
    set OutputName = "HADDrB_"
else
    echo "Error INVALID location!"
    exit
endif

set MainDir=`pwd`


set numFiles = 1
@ numFiles = ( $FileEnd - $FileStart ) / $FilesPerJob

set i = 0
while ($i <= $numFiles)
    echo "submitted"$i"/"$numFiles
    @ i = $i + 1
end

set XmlDir=./xml
if(! -e $XmlDir) exit
cd $XmlDir
set SubXml=sub.xml
if(-e $SubXml) rm $SubXml
touch $SubXml

set nFilePerJob=2 #40
set nFileTotal=all

# print xml file
echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${nFilePerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml
echo \<command\>$MainDir/run\.csh\</command\> >> $SubXml
# echo \<stdout URL=\"file:$MainDir/log/script_\$JOBINDEX\.out\" /\> >> $SubXml
echo \<stdout URL=\"file:/star/data01/pwg/svianping/log/script_\$JOBINDEX\.out\" /\> >> $SubXml
echo \<output fromScratch=\"root_\$JOBINDEX\.log\" toURL=\"file:/star/data01/pwg/svianping/log/\" /\> >> $SubXml
echo \<output fromScratch=\"output_\$JOBINDEX\.root\" toURL=\"file:/star/data01/pwg/svianping/output/\" /\> >> $SubXml
echo \<output fromScratch=\"KFParticleQA_\$JOBINDEX\.root\" toURL=\"file:/star/data01/pwg/svianping/output/\" /\> >> $SubXml
echo \</job\> >> $SubXml

star-submit $SubXml