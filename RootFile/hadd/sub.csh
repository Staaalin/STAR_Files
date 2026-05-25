#!/bin/csh

# For example:  ./sub.csh /star/data01/pwg/svianping/output/output_8 500

if ($#argv != 2) then
    echo "Usage: ./sub.csh <input_prefix> <files_per_job>"
    exit
endif

# set InputPrefix = $1
set InputPrefix = `echo $1 | sed 's/"//g'`
set FilesPerJob = $2

# DEBUG
echo "ARG1=[$1]"

# =========================
# count files
# =========================
set AllFiles = `ls ${InputPrefix}*.root | wc -l`
echo "Total files = $AllFiles"

set OutputDir = "/star/data01/pwg/svianping/hadd"
mkdir -p $OutputDir
cd $OutputDir

@ Start = 1
@ JobIndex = 1

set BadFileList = "BadFileList_${JobIndex}.list"
rm -f $BadFileList
touch $BadFileList

while ($Start <= $AllFiles)

    @ End = $Start + $FilesPerJob - 1
    if ($End > $AllFiles) then
        @ End = $AllFiles
    endif

    set SubXml = "hadd_\${JOBID}.xml"
    set FileList = "FileList_\${JOBID}.list"
    rm -f $SubXml
    touch $SubXml
    rm -f $FileList
    touch $FileList

    # =========================
    # XML header
    # =========================
    echo '<?xml version="1.0" encoding="utf-8" ?>' >> $SubXml
    echo '<job simulateSubmission="false" fileListSyntax="xrootd">' >> $SubXml

    # =========================
    # command
    # =========================
    echo \<shell\>singularity exec \-e \-B /direct \-B /star \-B /afs \-B /gpfs \-B /sdcc/lustre02 /cvmfs/star\.sdcc\.bnl\.gov/containers/rhic_sl7\.sif\</shell\> >> $SubXml # For a9
    echo \<command\> >> $SubXml
    # echo 'source setDEV2.csh' >> $SubXml

    echo 'hadd -f hadd_'\${JOBID}'.root @$FILELIST' >> $SubXml
    echo \</command\> >> $SubXml

    # =========================
    # SandBox（关键补充）
    # =========================

    # --- input ROOT files ---
    @ i = $Start

    while ($i <= $End)

        set f = "${InputPrefix}${i}.root"

        if ( -e $f ) then

            set size = `stat -c%s $f`

            if ( $size > 342 ) then
                echo $f >> $FileList
            else
                echo $f >> $BadFileList
            endif

        else
            echo $f >> $BadFileList
        endif

        @ i++

    end


    # =========================
    # output + logs
    # =========================
    echo \<input URL=\"filelist:$OutputDir/$FileList\" \/\> >> $SubXml
    echo \<stdout URL=\"file:${OutputDir}/hadd_\${JOBID}.log\" /\> >> $SubXml
    echo \<output fromScratch=\"hadd_\${JOBID}.root\" toURL=\"file:${OutputDir}/\" /\> >> $SubXml

    echo \</job\> >> $SubXml

    # =========================
    # submit
    # =========================
    star-submit-beta $SubXml

    echo "Submitted job $JobIndex : [$Start - $End]"

    @ Start = $End + 1
    @ JobIndex++

end

echo "Done."