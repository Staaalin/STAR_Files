#!/bin/csh


rm -rf /star/data01/pwg/svianping/QA/
mkdir /star/data01/pwg/svianping/QA/
mkdir /star/data01/pwg/svianping/QA/log/

set MainDir=`pwd`

set LeftBrackets = "\("
set RightBrackets = "\)"
set Quo = '\"'

set cen = 1
set opt_weight = 1



# 生成filelist
set input_file = "/star/u/svianping/STAR_Files/QA-Group/PSY/AuAu26p5_hpss_new.list"
set output_file = "modified_list.list"
set prefix = "root://xrdstar.rcf.bnl.gov:1095/"
set line_count = `wc -l < $input_file`


cd /star/data01/pwg/svianping/QA/
set i = 0
set j = 0
while ($i <= $numFiles)

    # set SubXml=sub.xml
    set SubXml=/star/data01/pwg/svianping/HADD/sub.xml
    if(-e $SubXml) rm $SubXml
    touch $SubXml

    # 读取第一行内容
    set first_line = `head -n 1 $input_file`

    # 添加前缀
    set modified_line = $prefix$first_line

    # 写入新文件
    echo $modified_line > $output_file

    set FILELIST = $output_file

    # print xml file
    echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
    echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${FilesPerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml

    echo \<command\> >> $SubXml
    echo rm $i\.log >> $SubXml
    # echo touch $i\.log >> $SubXml
    set ARM = " > "
    echo ll >> $SubXml
    echo echo \"000000000000000000000000000000000000000\" >> $SubXml
    echo root4star \-b RunAnalyzer_QA\.C$LeftBrackets$cen,$opt_weight,$Quo$FILELIST$Quo$RightBrackets
    echo ls  >> $SubXml
    echo \</command\> >> $SubXml

    echo \<SandBox installer=\"ZIP\"\> >> $SubXml
    echo \<Package name=\"ZIP\_File\_$i\"\> >> $SubXml
    # echo \<File\>file:/star/u/svianping/STAR\_Files/RootFile/HADDr\_xml\.C\</File\> >> $SubXml
    @ k = 0
    while ($k < $FilesPerJob)
        @ j = $FileStart + $i * $FilesPerJob + $k
        if ($j > $FileEnd) then
            break
        endif

        set FileName = $ObvInputName$j".root"
        if (-e $FileName) then
            echo \<File\>file:$FileName\</File\> >> $SubXml
        endif

        @ k = $k + 1
    end

    set MixEventPWD = "/star/u/svianping/STAR_Files/RootFile/mixevent/MixEvent.C"
    echo \<File\>file:$MixEventPWD\</File\> >> $SubXml

    echo \</Package\> >> $SubXml
    echo \</SandBox\> >> $SubXml
    echo \<stdout URL=\"file:/star/data01/pwg/svianping/HADD/log/script\_$i\.out\" /\> >> $SubXml
    echo \<output fromScratch=\"$i.log\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    set HC = "H_"
    set TC = "T_"
    echo \<output fromScratch=\"$OutputName$HC$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \<output fromScratch=\"$OutputName$TC$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \</job\> >> $SubXml

    star-submit $SubXml

    echo "submitted"$i"/"$numFiles
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/ZIP*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sch*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sub.xml
    # echo "Delate submit files"
    @ i = $i + 1
end
