#!/bin/csh


rm -rf /star/data01/pwg/svianping/QA/
mkdir /star/data01/pwg/svianping/QA/
mkdir /star/data01/pwg/svianping/QA/cen/
mkdir /star/data01/pwg/svianping/QA/log/
mkdir /star/data01/pwg/svianping/QA/xml/

set MainDir=`pwd`

set LeftBrackets = "\("
set RightBrackets = "\)"
set Quo = '\"'
set FilesPerJob = 1

set StoreL = " > Run_"
set StoreR = ".log "

set cen = 1
set opt_weight = 1

set OutputName = "cen1.v2.root"

# 生成filelist
set input_file = "/star/u/svianping/STAR_Files/QA-Group/PSY/AuAu26p5_hpss_new.list"
set line_count = `wc -l < $input_file`
set prefix = "root://xrdstar.rcf.bnl.gov:1095/"
set OutputURL = "/star/data01/pwg/svianping/QA/cen"
set OutputLogURL = "/star/data01/pwg/svianping/QA/log"

set numFiles = $line_count
cd /star/data01/pwg/svianping/QA/
set i = 0
set j = 0
echo "共有文件"$line_count"个"


while ($i <= $numFiles)


    # set SubXml=sub.xml
    set SubXml = /star/data01/pwg/svianping/QA/sub.xml
    if(-e $SubXml) rm $SubXml
    touch $SubXml
    echo file created

    set FILELIST = "modified_list_"$i".list"
    set output_file = "/star/data01/pwg/svianping/QA/"$FILELIST

    # 读取第一行内容
    # set first_line = `head -n $i $input_file`
    set first_line = `sed -n $i\p $input_file`

    # 添加前缀
    if ($i > 0) then
        echo prefix = $prefix
        echo first_line = $first_line
        echo modified_line = $prefix$first_line
    endif

    set modified_line = $prefix$first_line

    # 写入新文件
    echo $modified_line > $output_file

    echo FileList Created

    # print xml file
    echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
    echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${FilesPerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml

    echo \<command\> >> $SubXml
    
    echo setenv NODEBUG yes   >> $SubXml
    echo starver SL20d        >> $SubXml


    echo Environment Setted

    # echo touch $i\.log >> $SubXml
    set ARM = " > "
    echo ll >> $SubXml
    echo echo \"000000000000000000000000000000000000000\" >> $SubXml
    echo root4star \-b RunAnalyzer_QA\.C$LeftBrackets$cen,$opt_weight,$Quo$FILELIST$Quo$RightBrackets$StoreL$i$StoreR >> $SubXml
    echo ls  >> $SubXml
    echo mv $OutputName cen1_$i\.root >> $SubXml
    echo \</command\> >> $SubXml

    echo Command Setted

    echo \<SandBox installer=\"ZIP\"\> >> $SubXml
    echo \<Package name=\"ZIP\_File\_$i\"\> >> $SubXml
    # echo \<File\>file:/star/u/svianping/STAR\_Files/RootFile/HADDr\_xml\.C\</File\> >> $SubXml

    echo ZIP Start
    
    set StRootPWD = "/star/u/svianping/STAR_Files/QA-Group/PSY/StRoot"
    echo \<File\>file:$StRootPWD\</File\> >> $SubXml

    set GammaQAPWD = "/star/u/svianping/STAR_Files/QA-Group/PSY/Gamma_QA.C"
    echo \<File\>file:$GammaQAPWD\</File\> >> $SubXml

    set RunAnalyzerQAPWD = "/star/u/svianping/STAR_Files/QA-Group/PSY/RunAnalyzer_QA.C"
    echo \<File\>file:$RunAnalyzerQAPWD\</File\> >> $SubXml

    set ResolutionPWD = "/star/u/svianping/STAR_Files/QA-Group/PSY/Resolution_cen1.weight_112_QA_new.root"
    echo \<File\>file:$ResolutionPWD\</File\> >> $SubXml

    set FILELISTPWD = $output_file
    echo \<File\>file:$FILELISTPWD\</File\> >> $SubXml

    echo \<File\>file:$output_file\</File\> >> $SubXml

    echo \</Package\> >> $SubXml
    echo \</SandBox\> >> $SubXml

    echo SandBox Setted

    echo \<stdout URL=\"file:/star/data01/pwg/svianping/QA/log/script\_$i\.out\" /\> >> $SubXml
    echo \<output fromScratch=\"Run$i\.log\" toURL=\"file:$OutputLogURL\" /\> >> $SubXml
    echo \<output fromScratch=\"cen1\_$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \</job\> >> $SubXml

    cp $SubXml /star/data01/pwg/svianping/QA/xml/sub$i.xml
    star-submit $SubXml

    echo "submitted"$i"/"$numFiles

    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/ZIP*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sch*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sub.xml
    echo "Delate submit files"
    @ i = $i + 1


end


