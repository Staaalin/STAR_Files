#!/bin/csh

# 1. 保存 condor_q 输出
condor_q > temp.csh

# 2. 处理 temp.csh
awk '{
    # 跳过表头等非 job 行
    if ($1 !~ /^[0-9]+\./) next

    # 如果这一行里有字段等于大写 I，则跳过
    for (i=1; i<=NF; i++) {
        if ($i == "I") next
    }

    # 查找第一个以 .csh 结尾的字段
    for (i=1; i<=NF; i++) {
        if ($i ~ /\.csh$/) {
            jobid = $1
            script = $i
            # 把脚本后面所有参数拼起来
            cmd = $i
            for (j=i+1; j<=NF; j++) {
                cmd = cmd " " $j
            }
            print "condor_rm " jobid
            print "condor_submit " cmd " &"
            break
        }
    }
}' temp.csh > resubmit_jobs.csh

# 3. 给生成的脚本加执行权限
chmod +x resubmit_jobs.csh

echo "已生成 resubmit_jobs.csh，可以直接运行以批量 rm + submit"
