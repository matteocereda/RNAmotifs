# Configuration
# --------------------

# Before running this script, please adjust folders (paths) in the r.config and n.config files.

# SET THE ABSOLUTE PATH TO RNAMOTIFS FOLDER
RNAmotifs=$1 # SET THE ABSOLUTE PATH TO RNAMOTIFS FOLDER

cd $RNAmotifs

export PYTHONPATH=$RNAmotifs

# SET THE NAME OF THE RNAMOTIFS RESULTS FOLDER
result_dir=$2

# SET THE ABSOLUTE PATH TO SPLICING CHANGE PATH
splicing_change_infile=$RNAmotifs/input/$3

species=$4 # set species

organism=$5 # either "Mouse" or "Human"

p_emp=0.0005

p_fisher=0.1


# Run
# --------------------

N_config="${RNAmotifs}/input/${result_dir}.n.config"
R_config="${RNAmotifs}/input/${result_dir}.r.config"

# writing gMotifs configuration files

rm $N_config $R_config

if [ ! -f "$N_config" ] ; then
         # if not create the file
         touch "$N_config"
fi

echo "organism=$organism" >> $N_config
echo "splicing_change_file=$splicing_change_infile" >> $N_config
echo "tetramer_folder=$RNAmotifs/tetramers/$result_dir" >> $N_config
echo "results_folder=$RNAmotifs/results/$result_dir" >> $N_config

cp $N_config $R_config

echo "search=N" >> $N_config
echo "search=R" >> $R_config

cat $N_config

cat $R_config

# prepare non-overlapping regions from splicing file

if [ -e m3_light/regions/$result_dir ] ; then
	echo "Dir already present"
else
	mkdir m3_light/regions/$result_dir
fi

python m3_light/regions/prepare_regions.py $splicing_change_infile m3_light/regions/$result_dir/overlapping.tab

python m3_light/regions/merge_overlapping_regions.py m3_light/regions/$result_dir/overlapping.tab m3_light/regions/$result_dir/$result_dir.tab

rm m3_light/regions/$result_dir/overlapping.tab


# find tetramers in the genome

echo "import m3_light" > tmp.py
echo "m3_light.find_tetramers('${result_dir}', '${species}')" >> tmp.py
python tmp.py
rm tmp.py

# move results to tetramers folder

mkdir tetramers
mkdir tetramers/$result_dir
mkdir tetramers/$result_dir/r
mkdir tetramers/$result_dir/nr
cp m3_light/results/$result_dir/motifs.final/pth_0.5/A* tetramers/$result_dir/nr
cp m3_light/results/$result_dir/motifs.final/pth_0.5/C* tetramers/$result_dir/nr
cp m3_light/results/$result_dir/motifs.final/pth_0.5/T* tetramers/$result_dir/nr
cp m3_light/results/$result_dir/motifs.final/pth_0.5/G* tetramers/$result_dir/nr
cp m3_light/results/$result_dir/motifs.final/pth_0.5/W* tetramers/$result_dir/r
cp m3_light/results/$result_dir/motifs.final/pth_0.5/Y* tetramers/$result_dir/r
cp m3_light/results/$result_dir/motifs.final/pth_0.5/S* tetramers/$result_dir/r
cp m3_light/results/$result_dir/motifs.final/pth_0.5/R* tetramers/$result_dir/r


# counting tetramer occurences in the genome

# couting: count the number of exons with tetramer occurrencies in the region of interest r1,r2,r3
# tetramer: calculate the positional specific tetramer occurrences along the RNA splicing map

mkdir results
mkdir results/$result_dir
./gMotifs/build/tetramer -c $N_config
./gMotifs/build/counting -c $N_config
./gMotifs/build/tetramer -c $R_config
./gMotifs/build/counting -c $R_config

# bootstrap and tetramer selection

cd R

Rscript bootstrap-FDR.R ../results $result_dir 10000

Rscript selection_of_tetramers.R ../results $result_dir 10000 $p_fisher $p_emp
