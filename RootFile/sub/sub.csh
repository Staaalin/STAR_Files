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
    set OutputURL = "/star/data01/pwg/svianping/HADD/"
else if ($InputNameIndex == 2) then
    set ObvInputName = "/star/data01/pwg/svianping/HADD/HADD_"
    set ObvOutputName = "/star/data01/pwg/svianping/HADD/HADDrA_"
    set InputName = "HADD_"
    set OutputName = "HADDrA_"
    set OutputURL = "/star/data01/pwg/svianping/HADD/"
else if ($InputNameIndex == 3) then
    set ObvInputName = "/star/data01/pwg/svianping/HADD/HADDrA_"
    set ObvOutputName = "/star/u/svianping/STAR_Files/RootFile/HADDrB_"
    set InputName = "HADDrA_"
    set OutputName = "HADDrB_"
    set OutputURL = "/star/u/svianping/STAR_Files/RootFile/"
else
    echo "Error INVALID location!"
    exit
endif

set MainDir=`pwd`


set numFiles = 1
@ numFiles = ( $FileEnd - $FileStart ) / $FilesPerJob

cd /star/data01/pwg/svianping/HADD/
set i = 0
set j = 0
while ($i <= $numFiles)

    # set SubXml=sub.xml
    set SubXml=/star/data01/pwg/svianping/HADD/sub.xml
    if(-e $SubXml) rm $SubXml
    touch $SubXml

    # print xml file
    echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
    echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${FilesPerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml

    echo \<command\> >> $SubXml
    echo "source /star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh" >> $SubXml
    echo rm $i\.log >> $SubXml
    # echo touch $i\.log >> $SubXml
    set ARM = " > "
    echo ll >> $SubXml
    echo echo \"000000000000000000000000000000000000000\" >> $SubXml
    echo root4star -q -b \'HADDr_xml.C\(\"$InputName\",\"$OutputName\",$i,$FilesPerJob,$FileStart,$FileEnd\)\'$ARM$i".log" >> $SubXml
    # echo root4star -q -b \'HADDr_xml.C\(\"$InputName\",\"$OutputName\",$i,$FilesPerJob,$FileStart,$FileEnd\)\'" >> /star/u/svianping/STAR_Files/RootFile/sub.log" >> $SubXml
    echo \</command\> >> $SubXml

    echo \<SandBox installer=\"ZIP\"\> >> $SubXml
    echo \<Package name=\"ZIP\_File\_$i\"\> >> $SubXml
    echo \<File\>file:/star/u/svianping/STAR\_Files/RootFile/HADDr\_xml\.C\</File\> >> $SubXml
    set k = 0
    while ($k < $FilesPerJob)
        @ j = $FileStart + $i * $FilesPerJob + $k
        if ($j > $FileEnd) then
            break
        endif

        set FileName = $ObvInputName$j".root"
        if (-e $FileName) then
            echo \<File\>file:$FileName\</File\> >> $SubXml
        else
            break
        endif

        @ k = $k + 1
    end

    echo \</Package\> >> $SubXml
    echo \</SandBox\> >> $SubXml
    echo \<stdout URL=\"file:/star/data01/pwg/svianping/HADD/log/script\_$i\.out\" /\> >> $SubXml
    echo \<output fromScratch=\"$i.log\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \<output fromScratch=\"$OutputName$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \</job\> >> $SubXml

    star-submit $SubXml

    echo "submitted"$i"/"$numFiles
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/ZIP*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sch*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sub.xml
    # echo "Delate submit files"
    @ i = $i + 1
end
