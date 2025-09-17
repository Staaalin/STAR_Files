#!/bin/bash

# 先执行 condor_q > temp.log
rm temp.log
rm resubmit_jobs.csh
condor_q > temp.log

# 输出文件
outfile="resubmit_jobs.csh"
> "$outfile"   # 清空旧文件

# 遍历 temp.log
while read -r line; do
    # 判断是否以数字开头，且包含 " I " 和 ".csh"
    if [[ "$line" =~ ^[0-9]+\.[0-9]+ ]] && [[ "$line" == *" I "* ]] && [[ "$line" == *".csh"* ]]; then
        # 提取 jobid
        jobid=$(echo "$line" | awk '{print $1}')
        # 提取从第一个 .csh 开始的内容
        cmd=$(echo "$line" | sed -E 's/.* ([^ ]+\.csh .*)/\1/')

        echo "condor_rm $jobid" >> "$outfile"
        echo "condor_submit $cmd &" >> "$outfile"
    fi
done < temp.log

echo "已生成 $outfile"
