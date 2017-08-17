for file in $(ls -1 qq*.lhef)
do
    lastline=$(tail -1 $file)
    # echo $lastline
    if [ "$lastline" != "</LesHouchesEvents>" ]; then
        echo $file": fails"
    fi
done
