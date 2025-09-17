#!/bin/csh

# 1. 运行 condor_q 输出
condor_q > temp.csh

# 2. 用 awk 处理 temp.csh
awk '{
    # 如果这一行里有一个字段等于大写 R，则跳过
    for (i=1; i<=NF; i++) {
        if ($i == "R") next
    }

    # 否则正常处理
    if ($1 ~ /^[0-9.]+$/ && $NF ~ /\.csh$/) {
        jobid = $1
        script = $NF
        print "condor_rm " jobid
        print "condor_submit " script " &"
    }
}' temp.csh > resubmit_jobs.csh

# 3. 给生成的脚本加执行权限
chmod +x resubmit_jobs.csh

echo "已生成 resubmit_jobs.csh，可以直接运行以批量 rm + submit"
