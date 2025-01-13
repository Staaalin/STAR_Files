#!/bin/csh

# 定义目录路径
set dir = "/star/data01/pwg/svianping/log/"

# 循环从 script_0.out 到 script_10123.out
foreach i (`seq 0 10123`)
    # 生成文件名
    set filename = "${dir}script_${i}.out"
    
    # 检查文件是否存在
    if (! -e $filename) then
        # 文件不存在，创建新文件
        touch $filename
        echo "Created file: $filename"
    else
        # 文件已存在，跳过
        echo "File already exists: $filename"
    endif
end

# condor_q -hold -format '%s\n' ClusterId | xargs -I {} condor_release {}