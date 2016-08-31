mkdir hg19
cd hg19
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit 
../twoBitToFa hg19.2bit hg19.fasta -noMask
cd ../../
printf 'import m3_light\nm3_light.fastasplit("m3_light/genomes/hg19/hg19.fasta")' | python
mv *.string m3_light/genomes/hg19
rm m3_light/genomes/hg19/hg19.fasta
rm m3_light/genomes/hg19/hg19.2bit
