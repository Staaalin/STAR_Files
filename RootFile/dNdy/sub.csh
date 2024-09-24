#!/bin/csh

# set InputName = 
# set FilesPerJob = 400
# set FileStart = 1
# set FileEnd = 66389

echo "Particle PDG List:"
echo "+-2212    Proton"
echo "+-321     Kaon"
echo "+-211     Pion"
echo "  310     K0S"
echo "  333     Phi"
echo "+-3122    Lambda"
echo "+-3312    Xi"
echo "+-3334    Omega"

echo "Please enter particle A PDG:"
set A_PDG = "$<"
echo "Please enter particle B PDG:"
set B_PDG = "$<"


echo "Please enter which location:"
echo "SCHEME 1: /star/data01/pwg/svianping/output/output_*.root"
echo "SCHEME 2: /star/data01/pwg/svianping/HADD/HADD_T_*.root"
echo "SCHEME 3: /star/data01/pwg/svianping/HADD/HADDrA_*.root"
set InputNameIndex = "$<"

if ($InputNameIndex == 1) then

    cd /star/data01/pwg/svianping/output/
    set numFiles = `find . -maxdepth 1 -name "output_*.root" -type f | wc -l`

else if ($InputNameIndex == 2) then

    cd /star/data01/pwg/svianping/HADD/
    set numFiles = `find . -maxdepth 1 -name "HADD_T_*.root" -type f | wc -l`

endif

echo "一共有文件数：$numFiles"

echo "Please enter merge how much .root into ONE:"
set FilesPerJob = "$<"

echo "Set Start and End? 0:no , 1:yes"
set Mode = "$<"

if ($Mode == 1) then

    echo "Please enter from which file:"
    set FileStart = "$<"

    echo "Please enter to which file:"
    set FileEnd = "$<"

else if ($Mode == 0) then

    echo "Please enter MAX files scan:"
    set FileStart = 0
    set FileEnd = "$<"

endif


if ($InputNameIndex == 1) then
    set ObvInputName = "/star/data01/pwg/svianping/output/output_"
    set ObvOutputName = "/star/data01/pwg/svianping/HADD/HADD_"
    set InputName = "output_"
    set OutputName = "HADD_"
    set OutputURL = "/star/data01/pwg/svianping/HADD/"
    rm -rf /star/data01/pwg/svianping/HADD
    mkdir /star/data01/pwg/svianping/HADD
    mkdir /star/data01/pwg/svianping/HADD/log/
else if ($InputNameIndex == 2) then
    set ObvInputName = "/star/data01/pwg/svianping/HADD/HADD_T_"
    set ObvOutputName = "/star/data01/pwg/svianping/HADD/HADDr_"
    set InputName = "HADD_T_"
    set OutputName = "HADDr_"
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
    echo "source setDEV2.csh" >> $SubXml
    echo rm $i\.log >> $SubXml
    # echo touch $i\.log >> $SubXml
    set ARM = " > "
    echo ll >> $SubXml
    echo echo \"000000000000000000000000000000000000000\" >> $SubXml
    echo set midname = \"$InputName\" >> $SubXml
    echo set outmidname = \"$OutputName\" >> $SubXml
    set Jnum = 0
    @ Jnum = $FileStart + $i * $FilesPerJob
    echo set StartFileIndex = $Jnum >> $SubXml
    @ Jnum = $FileStart + ( $i + 1 ) * $FilesPerJob - 1
    echo set EndFileIndex = $Jnum >> $SubXml
    echo set OutputFileIndex = $i >> $SubXml
    echo set A_PDG = $A_PDG >> $SubXml
    echo set B_PDG = $B_PDG >> $SubXml
    set LeftBrackets = "\("
    set RightBrackets = "\)"
    set Quo = '\"'
    echo root4star \-b dNdy\.C$LeftBrackets$Quo\$midname$Quo,\$StartFileIndex,\$EndFileIndex,\$OutputFileIndex,$Quo\$outmidname$Quo,\$A_PDG,\$B_PDG$RightBrackets >> $SubXml
    # echo root4star -q -b \'HADDr_xml.C\(\"$InputName\",\"$OutputName\",$i,$FilesPerJob,$FileStart,$FileEnd\)\'$ARM$i".log" >> $SubXml
    echo ls  >> $SubXml
    echo \</command\> >> $SubXml

    echo \<SandBox installer=\"ZIP\"\> >> $SubXml
    echo \<Package name=\"ZIP\_File\_$i\"\> >> $SubXml
    # echo \<File\>file:/star/u/svianping/STAR\_Files/RootFile/HADDr\_xml\.C\</File\> >> $SubXml
    @ k = 0
    while ($k < $FilesPerJob)
        @ j = $FileStart + $i * $FilesPerJob + $k
        if ($j > $FileEnd) then
            break
        endif

        set FileName = $ObvInputName$j".root"
        if (-e $FileName) then
            echo \<File\>file:$FileName\</File\> >> $SubXml
        endif

        @ k = $k + 1
    end

    set MixEventPWD = "/star/u/svianping/STAR_Files/RootFile/dNdy/dNdy.C"
    echo \<File\>file:$MixEventPWD\</File\> >> $SubXml
    set SourceFilePWD = "/star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh"
    echo \<File\>file:$SourceFilePWD\</File\> >> $SubXml

    echo \</Package\> >> $SubXml
    echo \</SandBox\> >> $SubXml
    echo \<stdout URL=\"file:/star/data01/pwg/svianping/HADD/log/script\_$i\.out\" /\> >> $SubXml
    echo \<output fromScratch=\"$i.log\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    set HC = "H_"
    set TC = "T_"
    echo \<output fromScratch=\"$OutputName$HC$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \<output fromScratch=\"$OutputName$TC$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \</job\> >> $SubXml

    star-submit $SubXml

    echo "submitted"$i"/"$numFiles
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/ZIP*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sch*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sub.xml
    # echo "Delate submit files"
    @ i = $i + 1
end

echo This is $A_PDG - $B_PDG Corralation