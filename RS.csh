#!/bin/csh -f
# 禁用 history
set history = 0
unset savehist

# 1) 获取 condor_q 输出
rm temp.csh
rm resubmit_jobs.csh
condor_q > temp.csh

# 2) 输出文件
set outfile = resubmit_jobs.csh
echo "#!/bin/csh -f" > $outfile
echo "" >> $outfile

# 3) 逐行读取
foreach line ("`cat temp.csh`")
    # 按空格切分
    set words = ( $line )
    if ( $#words == 0 ) continue

    # 第一列必须是 jobid 格式：数字.数字
    if ("$words[1]" !~ [0-9]*"."[0-9]*) continue

    # 如果有字段等于 I，跳过
    set skip = 0
    foreach w ( $words )
        if ( "$w" == "I" ) then
            set skip = 1
        endif
    end
    if ( $skip == 1 ) continue

    # 找到第一个 .csh 脚本及后面的参数
    set idx = 0
    @ i = 1
    while ( $i <= $#words )
        if ( "$words[$i]" =~ "*.csh" ) then
            set idx = $i
            break
        endif
        @ i++
    end
    if ( $idx == 0 ) continue

    set cmd = ""
    @ j = $idx
    while ( $j <= $#words )
        set cmd = "$cmd $words[$j]"
        @ j++
    end

    echo "condor_rm $words[1]" >> $outfile
    echo "condor_submit$cmd &" >> $outfile
end

chmod +x $outfile
echo "已生成 $outfile"
