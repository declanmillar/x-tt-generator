#!/bin/bash

# check arguments
if (( $# != 2 ))
then
  echo "Usage:  add-lhefs.sh filename1 filename2"
  echo "Result: adds events in filename2 to filename1"
  exit 1
fi

# make newlines the only separator
IFS=$'\n'

# disable globbing
set -f

# declare count as an integer
declare -i count1
declare -i count2
declare -i count3

# count target events
count1=0
for i in $(cat < $1)
do
    if [ $i == "</event>" ]
    then
        count1+=1
    fi
done
echo "$1: $count1 events"
if [ $count1 = 10000 ]
then
    echo "target file already has 10,000 events; exiting"
    exit
fi

# count source events
count2=0
for i in $(cat < $2)
do
    if [ $i == "</event>" ]
    then
        count2+=1
    fi
done
echo "$2: $count2 events"
if [ $count2 = 10000 ]
then
    echo "source file has 10,000 events; exiting"
    exit
fi

count3=$(($count1 + $count2))
echo "$1 will have: $count3 events"
if (($count3 != 10000))
then
    echo "WARNING $1 will not have 10,000 events"
fi

if [ $count2 = 10000 ]
then
    echo "source file has 10,000 events; exiting"
    exit
fi

while true; do
    lastline=$(tail -1 $1 | head -1)
    # echo $lastline

    if [ "$lastline" == "</LesHouchesEvents>" ]
    then
        echo "deleting lhef closing tag in $1"
        sed '$ d' $1 > "$1.tmp"
        rm $1
        mv "$1.tmp" $1
        break
    elif [ "$lastline" == "</event>" ]
    then
        echo "target file ready; doing nothing"
        break
    else
        echo "event open; deleting last line"
        sed '$ d' $1 > "$1.tmp"
        rm $1
        mv "$1.tmp" $1
        continue
    fi
done

echo "copying events from $2 to $1"
copy=false
for i in $(cat $2)
do
    # echo $i
    if [ $i == "<event>" ]
    then
        copy=true
    fi
    if [ "$copy"==true ]
    then 
        echo "$i" >> $1
    fi
done

# count final events
count3=0
for i in $(cat < $1)
do
    if [ $i == "</event>" ]
    then
        count3+=1
    fi
done
echo "$1: $count3 events"

if [ $count3 = 10000 ]
then
    echo "copy success; deleting $2"
    rm $2
else
    echo "something went wrong: $2 not deleted"
fi