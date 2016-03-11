mkdir mm9
cd mm9
wget http://hgdownload-test.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.2bit -O mm9.2bit
# if not installed downaoload  http://hgdownload.cse.ucsc.edu/admin/exe/twoBitToFa 
../twoBitToFa mm9.2bit mm9.fasta -noMask
cd ../../
printf 'import m3_light\nm3_light.fastasplit("m3_light/genomes/mm9/mm9.fasta")' | python
mv *.string m3_light/genomes/mm9
rm m3_light/genomes/mm9/mm9.fasta
rm m3_light/genomes/mm9/mm9.2bit
