
# Parse arguments for input path, output path
while getopts 'x:c:' flag
do 
    case "${flag}" in 
        (x) PPIXpress_dir=${OPTARG};; # ../data/PPIXpress/INTACT_PPI/ResultFiles_renamed
        (c) PPICompare_dir=${OPTARG};; # ../data/PPICompare/INTACT_PPI
        (:) 
            case ${OPTARG} in 
                (x) exit 9999;;
                (c) exit 9999;;
            esac;;
    esac
done

mkdir -p $PPICompare_dir

get_ID (){
    ID=${1##*/};
    ID=${ID%%.*};
    echo $ID
}

for file1 in `ls $PPIXpress_dir`; 
do (
    for file2 in `ls $PPIXpress_dir`; 
    do (
        if [[ $file1 < $file2 ]]; then (
            ID1=`get_ID $file1`
            ID2=`get_ID $file2`

            echo $PPIXpress_dir/$file1 $PPIXpress_dir/$file2 $PPICompare_dir/${ID1}_${ID2}
            java -jar ../tools/PPICompare/PPICompare.jar -fdr=1 $PPIXpress_dir/$file1 $PPIXpress_dir/$file2 $PPICompare_dir/${ID1}_${ID2}
        ) fi
    ) 
    done
    ) 
done
echo "Finished"
