count_dir=data/Epispliced/htseq_count

get_ID (){
    ID=${1##*/};
    ID=${ID%%_*};
    echo $ID
}

for file1 in `ls $count_dir/*.txt`; 
do (
    for file2 in `ls $count_dir/*.txt`; 
    do (
        if [[ $file1 < $file2 ]]; then (
            ID1=`get_ID $file1`
            ID2=`get_ID $file2`
            Rscript code/normalize_count.R -f $count_dir -a $ID1 -b $ID2 -n 2
        ) fi
    ) 
    done
    ) 
done
echo "Finished"
