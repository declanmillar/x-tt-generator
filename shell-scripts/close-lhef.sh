#!/bin/bash

# make newlines the only separator
IFS = $'\n'

# disable globbing
set -f

# declare count as an integer
declare -i count

for file; do

    # count events
    count = 0

    for i in $(cat <$file); do
        # echo "tester: $i"
        if [ $i == "</event>" ]; then
            count += 1
        fi
    done

    echo "$file: $count events"

    # check file has at least one event
    if [ $count -lt 1 ]; then
        echo "skipping eventless file"
        continue
    fi

    while true; do
        lastline = $(tail -1 $file | head -1)
        # echo $lastline

        if [ "$lastline" == "</LesHouchesEvents>" ]; then
            echo "file already complete; doing nothing"
            break
        elif [ "$lastline" == "</event>" ]; then
            echo "</LesHouchesEvents>" >>$file
            echo "event closed; appended </LesHouchesEvents>"
        else
            echo "event open; deleting last line"
            sed '$ d' $file >"$file.tmp"
            rm $file
            mv "$file.tmp" $file
            continue
        fi
    done
done
