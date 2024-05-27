#!/bin/csh

rm -rf /star/data01/pwg/svianping/HADD/
mkdir /star/data01/pwg/svianping/HADD/

cd /star/data01/pwg/svianping/output/

# set numFiles = 'ls -l output_1000*.root | wc -l'
set numFiles = `find . -maxdepth 1 -name "output_*.root" -type f | wc -l`

echo "一共有文件数：$numFiles"

set i = 0
set j = 0
set k = 0
set n = 0
set start = 0
set end = 0
set step = 180

while ($n < 500)

    cd /star/data01/pwg/svianping/output/

    if (-e output_$i.root) then
        
        cp output_$i.root /star/data01/pwg/svianping/HADD/
        echo "成功复制 output_$i.root"
        @ n = 0

        @ j = $j + 1

    else
        @ n = $n + 1
    endif

    if ($j == $step) then
        cd /star/data01/pwg/svianping/HADD/
        hadd HADD_$k.root output_*.root
        # echo "成功创建 HADD_$k.root"
        rm -rf output_*.root
        @ k = $k + 1
        @ j = 0
    endif

    @ i = $i + 1


end


if (-e output_*.root) then
    hadd HADD_$k.root output_*.root
endif


set Energy = 62.0
set midname = "/star/data01/pwg/svianping/output/output_"
set StartFileIndex = 500
set EndFileIndex = 500
set OutputFileIndex = 0

root4star -b Subtract.C\(\"$midname\",$StartFileIndex,$EndFileIndex,$OutputFileIndex,$Energy\)