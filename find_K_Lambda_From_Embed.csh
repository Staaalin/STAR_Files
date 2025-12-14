#!/bin/csh

set BASE = /star/embed/embedding

foreach dir ($BASE/*)
    if (! -d $dir) continue

    set hasK = 0
    set hasLambda = 0

    # 查找包含 K 的子目录
    if (`ls -d $dir/*K* >& /dev/null; echo $?` == 0) then
        set hasK = 1
    endif

    # 查找包含 Lambda 的子目录
    if (`ls -d $dir/*Lambda* >& /dev/null; echo $?` == 0) then
        set hasLambda = 1
    endif

    # 同时满足
    if ($hasK && $hasLambda) then
        echo $dir
    endif
end
