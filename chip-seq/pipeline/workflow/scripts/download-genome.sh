#!/usr/bin/env bash
set -euo pipefail

cd $1

# tar.gz
wget -O genome.tar.gz https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
tar xvzf genome.tar.gz
rm genome.tar.gz

# single fasta file
rm -f mm10.fa
touch mm10.fa
for chr in *.fa
do
    cat $chr >> genome.fa
    rm -f $chr
done

