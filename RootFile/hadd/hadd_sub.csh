#!/bin/csh

if ($#argv != 2) then
    echo "Usage: ./hadd_sub.csh <input_prefix> <files_per_job>"
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

while ($Start <= $AllFiles)

    @ End = $Start + $FilesPerJob - 1
    if ($End > $AllFiles) then
        @ End = $AllFiles
    endif

    set SubXml = "hadd_${JobIndex}.xml"
    rm -f $SubXml
    touch $SubXml

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

    set FileList = ""
    @ i = $Start

    while ($i <= $End)
        set f = "${InputPrefix}${i}.root"
        if (-e $f) then
            set FileList = "$FileList $f"
        endif
        @ i++
    end

    # echo "hadd ${OutputDir}/hadd_${JobIndex}.root $FileList" >> $SubXml
    echo hadd hadd_${JobIndex}.root *.root >> $SubXml
    echo \</command\> >> $SubXml

    # =========================
    # SandBox（关键补充）
    # =========================
    echo \<SandBox installer="ZIP"\> >> $SubXml
    echo \<Package name=\"ZIP_File_${JobIndex}\"\> >> $SubXml

    # --- input ROOT files ---
    @ i = $Start
    while ($i <= $End)
        set f = "${InputPrefix}${i}.root"
        if (-e $f) then
            echo "<File>file:${f}</File>" >> $SubXml
        endif
        @ i++
    end

    # --- important: environment script ---
    # echo "<File>file:/star/u/svianping/STAR_Files/KFParticle4Lambda/setDEV2.csh</File>" >> $SubXml

    echo \</Package\> >> $SubXml
    echo \</SandBox\> >> $SubXml

    # =========================
    # output + logs
    # =========================
    echo \<stdout URL=\"file:${OutputDir}/hadd_${JobIndex}.log\" /\> >> $SubXml
    echo \<output fromScratch=\"hadd_${JobIndex}.root\" toURL=\"file:${OutputDir}/\" /\> >> $SubXml

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