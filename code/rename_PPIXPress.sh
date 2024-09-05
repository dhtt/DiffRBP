result_dir=../data/PPIXpress/INTACT_PPI
renamed_file_dir=$result_dir/ResultFiles_renamed
mkdir -p $renamed_file_dir

# Script to change file names
for file in $result_dir/ResultFiles/*; 
do (
        ID=${file##*/};
        ID=${ID%%_*};
        if [[ $ID != *'reference'* ]]; then (
            echo $file
            SAMPLE_NAME=`grep -B 1 ^${ID}_ppin.txt $result_dir/SampleSummary.txt | grep -v 'ppin'`
            mkdir -p $renamed_file_dir/${SAMPLE_NAME%%_*}
            mv $file $renamed_file_dir/${SAMPLE_NAME%%_*}/
            zip -r $renamed_file_dir/${SAMPLE_NAME%%_*}.zip $renamed_file_dir/${SAMPLE_NAME%%_*}
        ) 
        fi   
    )
done
echo "Finished"

