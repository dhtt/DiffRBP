
for file in ../data/Epispliced/raw_data/*.bam;
do (
    ID=${file%%.bam};
    ID=${ID##*/};
    htseq-count -f bam -t gene $file ../data/refgen/GCF_000001405.39_GRCh38.p13_genomic.gtf > $ID.txt
    )
done
wait 

# Script to change file names
for file in ../data/Epispliced/htseq_count/*; 
do (
    ID=${file%%.txt};
    ID=${ID##*/};
    FILENAME=$(grep $ID ../data/Epispliced/raw_data/rename.csv | awk '{print $3".txt"}'); 
    mv $file ../data/Epispliced/htseq_count/$FILENAME
    ) 
done
echo "Finished"