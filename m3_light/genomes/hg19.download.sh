mkdir hg19
cd hg19
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit -O hg19.2bit
../twoBitToFa hg19.2bit hg19.fasta -noMask
printf 'import m3_light\nm3_light.fastasplit("hg19.fasta")' | python
rm hg19.fasta
rm hg19.2bit
