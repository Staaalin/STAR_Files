#!/bin/csh

echo "Please enter how much datalists per job:"
set NDPJ = "$<"

rm -rf /star/data01/pwg/svianping/CD/
mkdir /star/data01/pwg/svianping/CD/
mkdir /star/data01/pwg/svianping/CD/cen/
mkdir /star/data01/pwg/svianping/CD/log/
mkdir /star/data01/pwg/svianping/CD/xml/

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
set OutputURL = "/star/data01/pwg/svianping/CD/cen"
set OutputLogURL = "/star/data01/pwg/svianping/CD/log"

set numFiles = $line_count
cd /star/data01/pwg/svianping/CD/
set i = 1
set j = 0
set k = 0
echo "共有文件"$line_count"个"

cp -r /star/u/svianping/STAR_Files/QA-Group/PSY/StRoot/StPicoEvent /star/u/svianping/STAR_Files/QA-Group/PSY/StPicoEvent

while ($i <= $numFiles)


    if ($j > $numFiles) then
    
        break

    endif

    # set SubXml=sub.xml
    set SubXml = /star/data01/pwg/svianping/CD/sub.xml
    if(-e $SubXml) rm $SubXml
    touch $SubXml
    echo file created

    set FILELIST = "modified_list_"$i".list"
    set output_file = "/star/data01/pwg/svianping/CD/"$FILELIST
    if(-e $output_file) rm $output_file
    touch $output_file

    @ k = 0
    while ($k < $NDPJ)

        if ($j > $numFiles) then

            break

        endif
        # 读取第一行内容
        # set first_line = `head -n $i $input_file`
        set first_line = `sed -n $j\p $input_file`

        # 添加前缀
        if ($i > 0) then
            echo prefix = $prefix
            echo first_line = $first_line
            echo modified_line = $prefix$first_line
        endif

        set modified_line = $prefix$first_line

        # 写入新文件
        echo $modified_line >> $output_file

        @ k = $k + 1

        @ j = $j + 1

    end

    echo FileList Created

    # print xml file
    echo \<\?xml version=\"1\.0\" encoding=\"utf-8\" \?\> >> $SubXml
    echo \<job simulateSubmission =\"false\" maxFilesPerProcess =\"${FilesPerJob}\" fileListSyntax=\"xrootd\"\> >> $SubXml
    echo \<shell\>singularity exec \-e \-B /direct \-B /star \-B /afs \-B /gpfs \-B /sdcc/lustre02 /cvmfs/star\.sdcc\.bnl\.gov/containers/rhic_sl7\.sif\</shell\> >> $SubXml # For a9

    echo \<command\> >> $SubXml
    
    echo setenv NODEBUG yes   >> $SubXml
    echo starver SL24c        >> $SubXml
    # echo starver SL20d        >> $SubXml


    echo Environment Setted

    # echo touch $i\.log >> $SubXml
    set ARM = " > "
    echo mkdir StRoot >> $SubXml
    echo mv StPicoEvent StRoot/ >> $SubXml
    echo ll >> $SubXml
    echo echo \"000000000000000000000000000000000000000\" >> $SubXml
    echo root4star \-b \-q RunCD\.C$LeftBrackets$Quo$FILELIST$Quo$RightBrackets$StoreL$i$StoreR >> $SubXml
    # echo root4star \-b \-q CD\.C$LeftBrackets$Quo$FILELIST$Quo$RightBrackets$StoreL$i$StoreR >> $SubXml
    echo ls  >> $SubXml
    echo mv $OutputName cen1_$i\.root >> $SubXml
    echo \</command\> >> $SubXml

    echo Command Setted

    echo \<SandBox installer=\"ZIP\"\> >> $SubXml
    echo \<Package name=\"ZIP\_File\_$i\"\> >> $SubXml
    # echo \<File\>file:/star/u/svianping/STAR\_Files/RootFile/HADDr\_xml\.C\</File\> >> $SubXml

    echo ZIP Start
    
    set StRootPWD = "/star/u/svianping/STAR_Files/QA-Group/PSY/StPicoEvent"
    echo \<File\>file:$StRootPWD\</File\> >> $SubXml

    set GammaQAPWD = "/star/u/svianping/STAR_Files/QA-Group/PSY/CD.C"
    echo \<File\>file:$GammaQAPWD\</File\> >> $SubXml

    set RunAnalyzerQAPWD = "/star/u/svianping/STAR_Files/QA-Group/PSY/RunCD.C"
    echo \<File\>file:$RunAnalyzerQAPWD\</File\> >> $SubXml

    set FILELISTPWD = $output_file
    echo \<File\>file:$FILELISTPWD\</File\> >> $SubXml

    echo \<File\>file:$output_file\</File\> >> $SubXml

    echo \</Package\> >> $SubXml
    echo \</SandBox\> >> $SubXml

    echo SandBox Setted

    echo \<stdout URL=\"file:/star/data01/pwg/svianping/CD/log/script\_$i\.out\" /\> >> $SubXml
    echo \<output fromScratch=\"Run\_$i\.log\" toURL=\"file:$OutputLogURL\" /\> >> $SubXml
    echo \<output fromScratch=\"cen1\_$i\.root\" toURL=\"file:$OutputURL\" /\> >> $SubXml
    echo \</job\> >> $SubXml

    cp $SubXml /star/data01/pwg/svianping/CD/xml/sub$i.xml
    star-submit-beta $SubXml

    echo "submitted"$i"/"$numFiles

    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/ZIP*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sch*
    # rm -rf /star/u/svianping/STAR_Files/RootFile/sub/sub.xml
    echo "Delate submit files"
    @ i = $i + 1


end

rm -rf /star/u/svianping/STAR_Files/QA-Group/PSY/StPicoEvent

