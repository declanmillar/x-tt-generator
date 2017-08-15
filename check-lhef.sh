for file in $(ls -1 gg-tt-bbmumuvv_SM_13TeV_CT14LL_000.lhef*)
do
    lastline=$(tail -1 $file)
    # echo $lastline
    if [ "$lastline" != "</LesHouchesEvents>" ]; then
        echo $file, ": fails"
    fi
done
