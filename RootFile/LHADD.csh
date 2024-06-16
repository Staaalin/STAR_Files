#!/bin/csh
echo "准备执行"
sleep 5
echo "执行"

rm condor.log
condor_q > condor.log
# condor_rm svianping
# sleep 10

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
set start = 65204
set end = 0
set step = 100

@ i = $i + $start

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

    echo $j

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

    @ k = $k + 1

    hadd HADD_$k.root output_*.root

endif

# while ($i <= $numFiles)

#     cp output_$i.root /star/data01/pwg/svianping/HADD/
#     echo "成功复制 output_$i.root"

#     if ($j == $step) then
#         cd /star/data01/pwg/svianping/HADD/
#         hadd HADD_$k.root output_*.root
#         echo "成功创建 HADD_$k.root"
#         rm -rf output_*.root
#         @ j = 0
#     endif

#     @ i = $i + 1

#     @ j = $j + 1

# end

# @ numFiles = `find . -maxdepth 1 -name "HADD_*.root" -type f | wc -l`
# @ i = 0
# @ j = 0
# @ k = 0
# @ start = 0
# @ end = 0
# @ step = 50
# while ($i <= $numFiles)

#     cp HADD_$i.root /star/data01/pwg/svianping/output/
#     echo "成功复制 HADD_$i.root"

#     if ($j == $step) then
#         cd /star/data01/pwg/svianping/output/
#         hadd HADDC_$k.root HADD_*.root
#         echo "成功创建 HADDC_$k.root"
#         rm -rf HADD_*.root
#         @ j = 0
#     endif

#     @ i = $i + 1

#     @ j = $j + 1

# end