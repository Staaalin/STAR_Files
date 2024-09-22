#!/bin/csh

echo "Please enter output file index:"
set i = "$<"

rm -rf /star/data01/pwg/svianping/Random/
mkdir /star/data01/pwg/svianping/Random/

cd /star/data01/pwg/svianping/Random/

set FilesPerJob = 1

set OutputName = "Random_"
set OutputURL = /star/data01/pwg/svianping/Random/


# set SubXml=sub.xml
set SubXml=/star/data01/pwg/svianping/Random/sub.xml
if(-e $SubXml) rm $SubXml
touch $SubXml

# print xml file
echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${FilesPerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml

echo \<command\> >> $SubXml
echo "source setDEV2.csh" >> $SubXml
echo rm $i\.log >> $SubXml
# echo touch $i\.log >> $SubXml
set ARM = " > "
echo ll >> $SubXml
echo echo \"000000000000000000000000000000000000000\" >> $SubXml
set LeftBrackets = "\("
set RightBrackets = "\)"
set Quo = '\"'
echo root4star \-b RandomGenerator\.C$LeftBrackets$i$RightBrackets >> $SubXml
# echo root4star -q -b \'HADDr_xml.C\(\"$InputName\",\"$OutputName\",$i,$FilesPerJob,$FileStart,$FileEnd\)\'$ARM$i".log" >> $SubXml
echo ls  >> $SubXml
echo \</command\> >> $SubXml

echo \<SandBox installer=\"ZIP\"\> >> $SubXml
echo \<Package name=\"ZIP\_File\_$i\"\> >> $SubXml

set MixEventPWD = "/star/u/svianping/STAR_Files/RootFile/mixevent/RandomGenerator.C"
echo \<File\>file:$MixEventPWD\</File\> >> $SubXml
set SourceFilePWD = "/star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh"
echo \<File\>file:$SourceFilePWD\</File\> >> $SubXml

echo \</Package\> >> $SubXml
echo \</SandBox\> >> $SubXml
echo \<stdout URL=\"file:/star/data01/pwg/svianping/HADD/log/script\_$i\.out\" /\> >> $SubXml
echo \<output fromScratch=\"$i.log\" toURL=\"file:$OutputURL\" /\> >> $SubXml
set HC = "H_"
echo \<output fromScratch=\"$OutputName$HC$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
echo \</job\> >> $SubXml

star-submit $SubXml

echo "submitted"$i"/"$numFiles
# rm -rf /star/u/svianping/STAR_Files/RootFile/sub/ZIP*
# rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sch*
# rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sub.xml
# echo "Delate submit files"
