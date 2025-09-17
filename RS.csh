#!/bin/csh -f
# resubmit.csh  -- 生成 resubmit_jobs.csh（不使用 awk，不触发 Event not found）

# 1) 获取 condor_q 输出
rm temp.csh
rm resubmit_jobs.csh
condor_q > temp.csh

# 2) 输出文件初始化
set outfile = resubmit_jobs.csh
echo "#!/bin/csh -f" > $outfile
echo "" >> $outfile

# 3) 逐行读取并处理 temp.csh
while ( 1 )
    set line = "$<"
    if ( $status != 0 ) then
        break
    endif

    # 拆成字段数组
    set words = ( $line )
    if ( $#words == 0 ) continue

    # 只处理第一字段是 jobid（数字.数字）的行
    echo "$words[1]" | /bin/grep -E -q '^[0-9]+\.[0-9]+$'
    if ( $status != 0 ) continue

    # 如果任一字段是 I，则跳过
    set skip = 0
    @ i = 1
    while ( $i <= $#words )
        if ( "$words[$i]" == "I" ) then
            set skip = 1
            break
        endif
        @ i++
    end
    if ( $skip == 1 ) continue

    # 找到第一个以 .csh 结尾的字段
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

    # 拼接从 idx 到行尾的字段
    set cmd = "$words[$idx]"
    @ j = $idx + 1
    while ( $j <= $#words )
        set cmd = "$cmd $words[$j]"
        @ j++
    end

    # 写入 resubmit_jobs.csh
    echo "condor_rm $words[1]" >> $outfile
    echo "condor_submit $cmd &" >> $outfile
end < temp.csh

# 4) 可执行权限
chmod +x $outfile

echo "已生成 $outfile"
