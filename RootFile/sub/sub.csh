#!/bin/csh

set MainDir=`pwd`

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