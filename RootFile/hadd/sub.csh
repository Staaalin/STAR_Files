#!/bin/csh

if ($#argv < 2) then
    echo "Usage: hadd_sub /path/to/files/*.root FilesPerJob"
    exit 1
endif

set Pattern = $1
set FilesPerJob = $2

# 输出目录
set OutDir = "/star/data01/pwg/svianping/hadd/"
mkdir -p $OutDir
mkdir -p $OutDir/log

# 获取文件列表
set FileList = (`ls $Pattern`)
set TotalFiles = $#FileList

echo "Total files: $TotalFiles"
echo "Files per job: $FilesPerJob"

@ numJobs = ( $TotalFiles + $FilesPerJob - 1 ) / $FilesPerJob

echo "Total jobs: $numJobs"

@ i = 0
while ($i < $numJobs)

    set SubXml = "$OutDir/sub_$i.xml"
    if (-e $SubXml) rm $SubXml
    touch $SubXml

    # ---------- XML header ----------
    echo '<?xml version="1.0" encoding="utf-8" ?>' >> $SubXml
    echo '<job simulateSubmission="false" maxFilesPerProcess="'$FilesPerJob'" fileListSyntax="xrootd">' >> $SubXml

    # ---------- command ----------
    echo \<shell\>singularity exec \-e \-B /direct \-B /star \-B /afs \-B /gpfs \-B /sdcc/lustre02 /cvmfs/star\.sdcc\.bnl\.gov/containers/rhic_sl7\.sif\</shell\> >> $SubXml # For a9
    echo '<command>' >> $SubXml

    echo "echo START JOB $i" >> $SubXml

    set OutputFile = "hadd_$i.root"

    # 构建 hadd 命令
    set cmd = "hadd $OutputFile *.root"

    @ start = $i * $FilesPerJob
    @ end = $start + $FilesPerJob - 1

    if ($end >= $TotalFiles) then
        @ end = $TotalFiles - 1
    endif

    @ j = $start
    while ($j <= $end)
        set fname = $FileList[$j+1]
        set base = `basename $fname`
        set cmd = "$cmd $base"
        @ j++
    end

    echo $cmd >> $SubXml
    echo "ls" >> $SubXml

    echo '</command>' >> $SubXml

    # ---------- Resource ----------
    echo '<ResourceUsage>' >> $SubXml
    echo '<Priority>75</Priority>' >> $SubXml
    echo '</ResourceUsage>' >> $SubXml

    # ---------- Sandbox ----------
    echo '<SandBox installer="ZIP">' >> $SubXml
    echo '<Package name="ZIP_File_'$i'">' >> $SubXml

    @ j = $start
    while ($j <= $end)
        set fname = $FileList[$j+1]
        echo "<File>file:$fname</File>" >> $SubXml
        @ j++
    end

    echo '</Package>' >> $SubXml
    echo '</SandBox>' >> $SubXml

    # ---------- 输出 ----------
    echo '<stdout URL="file:'$OutDir'/log/job_'$i'.out" />' >> $SubXml
    echo '<output fromScratch="'$OutputFile'" toURL="file:'$OutDir'/" />' >> $SubXml

    echo '</job>' >> $SubXml

    # ---------- 提交 ----------
    star-submit-beta $SubXml

    echo "Submitted job $i / $numJobs"

    @ i++
end

echo "All jobs submitted."