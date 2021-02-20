#!/bin/bash

# make newlines the only separator
IFS = $'\n'

# disable globbing
set -f

# declare count as an integer
declare -i count

for file; do
    count = 0

    for i in $(cat <$file); do
        # echo "i: $i"
        if [ $i == "</event>" ]; then
            count += 1
        fi
    done

    echo "$file: $count events"
done
