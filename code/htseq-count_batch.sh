
for i in `ls ../data/Epispliced/raw_data/*.bam`
do 
    echo $i
done
# htseq-count -f bam -t gene test.bam ../data/refgen/GCF_000001405.39_GRCh38.p13_genomic.gtf > test.txt