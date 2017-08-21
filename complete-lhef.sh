#!/bin/bash

for file
do
    while true; do
        lastline=$(tail -1 $file | head -1)
        # echo $lastline

        if [ "$lastline" == "</LesHouchesEvents>" ]
        then
            echo "file already complete; doing nothing"
            break
        elif [ "$lastline" == "</event>" ]
        then
            echo "</LesHouchesEvents>" >> $file
            echo "event closed; appended </LesHouchesEvents>"
        else
            echo "event open; deleting last line"
            sed '$ d' $file > "$file.tmp"
            rm $file
            mv "$file.tmp" $file
            continue
        fi
    done
done
