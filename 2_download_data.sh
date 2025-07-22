#/usr/bin/env bash

# This script uses fasterq-dump to download SRA accessions listed in the input accession file

cpus=24

#accessions = $(ls ../accessions/)
#$accessions
country_data=../accessions/accession_USA3.txt


#cat $accessions | xargs fasterq-dump -e $cpus --outdir $output_dir
#for file in `ls ../accessions/accession_china.txt`
#do
f=$(basename $country_data .txt)
`mkdir -pv ../data/$f`
for accession in $(awk '!/run_accession/ {print $1}' $country_data)
 do
  echo downloading $accession

 fasterq-dump $accession -e $cpus -O ../data/$f
  echo download of the sample ended
 done
#done
