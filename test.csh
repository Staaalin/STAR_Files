#!/bin/csh
echo "开始执行"
sleep 7200
echo "执行"

condor_rm svianping

set ObvInputName = "/star/data01/pwg/svianping/output/output_"
set ObvOutputName = "/star/data01/pwg/svianping/HADD/HADD_"
set InputName = "output_"
set OutputName = "HADD_"
set OutputURL = "/star/data01/pwg/svianping/HADD/"

set FilesPerJob = 171
set FileStart = 0
set FileEnd = 68191

set MainDir=`pwd`


set numFiles = 1
@ numFiles = ( $FileEnd - $FileStart ) / $FilesPerJob

rm -rf /star/data01/pwg/svianping/HADD/*
mkdir /star/data01/pwg/svianping/HADD/log/
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
