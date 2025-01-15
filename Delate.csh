#!/bin/csh

rm -rf /star/data01/pwg/svianping/output/
rm -rf /star/data01/pwg/svianping/log/
rm -rf /star/data01/pwg/svianping/JobID/

mkdir /star/data01/pwg/svianping/output/
mkdir /star/data01/pwg/svianping/log/
mkdir /star/data01/pwg/svianping/JobID/ 

# 定义目录路径
set dir = "/star/data01/pwg/svianping/log/"

set j = 0
# 循环从 script_0.out 到 script_10123.out
foreach i (`seq 0 10123`)
    # 生成文件名
    set filename = "${dir}script_${i}.out"
    
    touch $filename

    @ j = $j + 1
    if ($j > 500) then
        echo "Done!"
        @ j = 0
    endif
end

echo "Finish Create Script File"
