#!/bin/csh

cd /star/data01/pwg/svianping/Subtract/

set i = 0
set j = 0
while (2 > 1)
    if (-e Cor_$i.root) then
        rm -rf /star/data01/pwg/svianping/Subtract/$i/
    endif
    if ($i < 20000) then
        @ i = $i + 1
    else
        sleep 60
        @ i = 0
        @ j = $j + 1
    endif
    if ($j > 720) then
        break
    endif
end