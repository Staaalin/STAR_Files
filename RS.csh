#!/bin/csh -f
# resubmit.csh  -- 生成 resubmit_jobs.csh（不使用 awk）

# 1) 获取 condor_q 输出
rm temp.csh
rm resubmit_jobs.csh
condor_q > temp.csh

# 2) 输出文件初始化
set outfile = resubmit_jobs.csh
echo "#!/bin/csh -f" >! $outfile
echo "" >>! $outfile

# 3) 逐行读取并处理 temp.csh
while ( 1 )
    set line = "$<"
    if ( $status != 0 ) then
        break
    endif

    # 拆成字段数组（按空白分割）
    set words = ( $line )
    if ( $#words == 0 ) then
        continue
    endif

    # 只处理以数字开头的第一字段（job 行），否则跳过
    echo "$words[1]" | /bin/grep -E -q '^[0-9]'
    if ( $status != 0 ) then
        continue
    endif

    # 如果任一字段精确等于大写 I，则跳过该行
    set skip = 0
    @ i = 1
    while ( $i <= $#words )
        if ( "$words[$i]" == "I" ) then
            set skip = 1
            break
        endif
        @ i = $i + 1
    end
    if ( $skip == 1 ) then
        continue
    endif

    # 寻找第一个以 .csh 结尾的字段（脚本名），并拼接它到行尾作为 submit 的参数
    set idx = 0
    @ i = 1
    while ( $i <= $#words )
        if ( "$words[$i]" =~ "*.csh" ) then
            set idx = $i
            break
        endif
        @ i = $i + 1
    end

    # 若找不到 .csh 字段则跳过
    if ( $idx == 0 ) then
        continue
    endif

    # 拼接从 idx 到行尾的所有字段为 cmd
    set cmd = "$words[$idx]"
    @ j = $idx + 1
    while ( $j <= $#words )
        set cmd = "$cmd $words[$j]"
        @ j = $j + 1
    end

    # 写入输出脚本
    echo "condor_rm $words[1]" >>! $outfile
    echo "condor_submit $cmd &" >>! $outfile

end < temp.csh

# 4) 赋予可执行权限
chmod +x $outfile

echo "已生成 $outfile 。运行 ./resubmit_jobs.csh 可执行批量 rm + submit（请先检查确认）"
