if [ $# != 2 ] echo "usage: check-lhef.sh <wildcard>"

for file in $(ls -1 $1)
do
    lastline=$(tail -1 $file)
    # echo $lastline
    if [ "$lastline" != "</LesHouchesEvents>" ]; then
        echo $file": fails"
    fi
done
