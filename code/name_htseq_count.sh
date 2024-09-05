# Script to change file names
for file in ../data/Epispliced/htseq_count/*; 
do (
    ID=${file%%.txt};
    ID=${ID##*/};
    FILENAME=$(grep $ID ../data/Epispliced/raw_data/rename.csv | awk '{print $3".txt"}'); echo $file ../data/Epispliced/htseq_count/$FILENAME
    ) 
done;
