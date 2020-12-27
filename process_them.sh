#!/usr/bin/env bash

GTF_FILE="${1%.*}"

# grep hashtag lines

mkdir temp_files;
mkdir results;
grep "#" $1 > ${GTF_FILE}_hashtag_lines
# split input gtf into fw and rev
awk -F'\t' '$7=="+" {print $0}' $1 > ${GTF_FILE}_fw.gtf
awk -F'\t' '$7=="-" {print $0}' $1 > ${GTF_FILE}_rev.gtf

for threshold in "500" "1000" "1500" "2000" "2500" "3000" "3500" "4000" "4500" "5000"; do
	python3 Extend_by_threshold_fw.py --input ${GTF_FILE}_fw.gtf --threshold $threshold > temp_files/${GTF_FILE}_fw.gtf ;
	python3 Extend_by_threshold_rev.py --input ${GTF_FILE}_rev.gtf --threshold $threshold > temp_files/${GTF_FILE}_rev.gtf ;
	cat temp_files/${GTF_FILE}_fw.gtf > temp_files/unsorted ;
	cat temp_files/${GTF_FILE}_rev.gtf >> temp_files/unsorted ;
	sort -k1,1V -k4,4n -k5,5rn temp_files/unsorted > temp_files/sorted ;
	cat ${GTF_FILE}_hashtag_lines > results/${GTF_FILE}_ext_by_${threshold}.gtf ;
	cat temp_files/sorted >> results/${GTF_FILE}_ext_by_${threshold}.gtf ;
	rm -rf temp_files/* ;
done

rm ${GTF_FILE}_hashtag_lines ;
rm ${GTF_FILE}_fw.gtf ;
rm ${GTF_FILE}_rev.gtf ;

