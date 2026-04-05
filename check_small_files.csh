#!/bin/csh

# 编写一个csh脚本，使其遍历/star/data01/pwg/svianping/output/路径下文件，
# 删除/显示所有大小小于9000的文件，这是两个模式：删除/显示，显示的话只显示头
# 五个以做检查，因为扫描数量很大不能全显示.

# =========================
# 参数检查
# =========================
if ( $#argv != 1 ) then
    echo "Usage: $0 [show|delete]"
    exit 1
endif

set mode = $1
set dir = /star/data01/pwg/svianping/output

# =========================
# 查找小文件（<9000 bytes）
# =========================
set files = `find $dir -type f -size -9000c`

if ( "$files" == "" ) then
    echo "No files smaller than 9000 bytes found."
    exit 0
endif

# =========================
# 模式：show
# =========================
if ( "$mode" == "show" ) then
    echo "Showing first 5 small files (<9000 bytes):"

    set count = 0
    foreach f ( $files )
        echo $f
        @ count++
        if ( $count >= 5 ) break
    end

# =========================
# 模式：delete
# =========================
else if ( "$mode" == "delete" ) then
    echo "Deleting all files smaller than 9000 bytes..."

    foreach f ( $files )
        rm -f "$f"
    end

    echo "Done."

else
    echo "Invalid mode: $mode"
    echo "Usage: $0 [show|delete]"
    exit 1
endif