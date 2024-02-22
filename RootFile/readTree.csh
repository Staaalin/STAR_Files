#!/bin/csh

set MainDir=`pwd`

set XmlDir=./xml
if(! -e $XmlDir) exit
cd $XmlDir
set SubXml=sub.xml
if(-e $SubXml) rm $SubXml
touch $SubXml

set nFilePerJob=3 #40
set nFileTotal=all

# print xml file
echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${nFilePerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml
echo \<command\>/star/u/svianping/STAR_Files/RootFile/readTree.C\</command\> >> $SubXml
echo \<stdout URL=\"file:$MainDir/log/script_readTree\.out\" /\> >> $SubXml
echo \</job\> >> $SubXml

star-submit $SubXml
