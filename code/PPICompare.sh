PPIXpress_dir=../data/PPIXpress/INTACT_PPI
PPICompare_dir=../data/PPICompare/INTACT_PPI
renamed_file_dir=$PPIXpress_dir/ResultFiles_renamed

mkdir -p $PPICompare_dir

get_ID (){
    ID=${1##*/};
    ID=${ID%%.*};
    echo $ID
}

for file1 in `ls $renamed_file_dir`; 
do (
    for file2 in `ls $renamed_file_dir`; 
    do (
        if [[ $file1 < $file2 ]]; then (
            ID1=`get_ID $file1`
            ID2=`get_ID $file2`
            java -jar ../../tools/PPICompare/PPICompare.jar $renamed_file_dir/$file1 $renamed_file_dir/$file2 $PPICompare_dir/${ID1}_${ID2}
        ) fi
    ) &
    done
    ) 
done
echo "Finished"
