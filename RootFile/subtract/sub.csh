#!/bin/csh

rm -rf /star/data01/pwg/svianping/Subtract/
mkdir /star/data01/pwg/svianping/Subtract/
mkdir /star/data01/pwg/svianping/Subtract/log/

cd /star/data01/pwg/svianping/output/

# set numFiles = 'ls -l output_1000*.root | wc -l'
# set numFiles = `find . -maxdepth 1 -name "output_*.root" -type f | wc -l`

# echo "一共有文件数：$numFiles"


set Energy = 62.0
set midname = "/star/data01/pwg/svianping/output/output_"

set i = 0
set j = 0
set k = 0
set n = 0
set step = 12
set StartFileIndex = 0
set EndFileIndex = 0
set OutputFileIndex = 0
set SubXml=sub.xml

cd /star/data01/pwg/svianping/output/
while ($n < 500)
    echo "Here good"
    if (($j == 0) && !(-e $SubXml)) then
        # 写入xml
        mkdir /star/data01/pwg/svianping/Subtract/$k/
        cd /star/data01/pwg/svianping/Subtract/$k/
        set SubXml = /star/data01/pwg/svianping/Subtract/$k/sub.xml
        touch $SubXml
        echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
        echo \<SandBox installer=\"ZIP\"\> >> $SubXml
        echo \<Package name=\"ZIP\_File\_$k\"\> >> $SubXml
        echo \<File\>file:/star/u/svianping/STAR\_Files/RootFile/subtract/Subtract\.C\</File\> >> $SubXml
    endif

    if (-e output_$i.root) then
        echo \<File\>file:/star/data01/pwg/svianping/output/output\_$i\.root\</File\> >> $SubXml
        @ n = 0

        @ j = $j + 1

    else
        @ n = $n + 1
    endif

    if (($j == $step) || (($n == 500) && ($j > 0))) then
        @ EndFileIndex = $i
        echo \</Package\> >> $SubXml
        echo \</SandBox\> >> $SubXml
        echo \<command\> >> $SubXml
        echo "source /star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh" >> $SubXml
        echo rm $i\.log >> $SubXml
        set ARM = " > "
        echo ll >> $SubXml
        echo echo \"000000000000000000000000000000000000000\" >> $SubXml
        echo root4star -q -b \'Subtract.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,$Energy\)\'$ARM$k".log" >> $SubXml
        echo \</command\> >> $SubXml
        echo \<stdout URL=\"file:/star/data01/pwg/svianping/Subtract/log/script\_$k\.out\" /\> >> $SubXml
        echo \<output fromScratch=\"$k.log\" toURL=\"file:/star/data01/pwg/svianping/Subtract/log/root\_$k\.log\" /\> >> $SubXml
        echo \<output fromScratch=\"Cor\_$k\.root\" toURL=\"file:/star/data01/pwg/svianping/Subtract/Cor\_$k\.root\" /\> >> $SubXml
        echo \</job\> >> $SubXml

        star-submit $SubXml
        
        @ StartFileIndex = $i + 1
        @ j = 0

        @ k = $k + 1

    endif

    @ i = $i + 1


end


