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
set start = 0
set end = 0
set step = 50
while ($i <= $numFiles)

    cp output_$i.root /star/data01/pwg/svianping/HADD/
    echo "成功复制 output_$i.root"

    if ($j == $step) then
        cd /star/data01/pwg/svianping/HADD/
        hadd HADD_$k.root output_*.root
        echo "成功创建 HADD_$k.root"
        rm -rf output_*.root
        @ j = 0
        break
    endif

    @ i = $i + 1

    @ j = $j + 1

end
