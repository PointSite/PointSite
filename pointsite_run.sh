#!/bin/bash

if [ $# -lt 4 ]
then
	echo "Usage: ./pointsite_run.sh <data> <indir> <outdir> <home> "
	echo "[note]: <data> shall be list containing the input PDB names, such as '1pazA'. "
	echo "        <indir> shall be the input directory that contains these input PDB files. "
	echo "        <outdir> shall be the output directory for the result. "
	echo "        <home> shall be the home directory of the program. "
	exit
fi

#----- input arguments ---------#
data=$1
indir=$2
outdir=$3
home=$4

#----- create temporary folder ----#
mkdir -p $outdir
outdir=`readlink -f $outdir`
fulnam=`basename $data`
relnam=${fulnam%.*}
tmpdir=$outdir/tmp_${relnam}_xyz
mkdir -p $tmpdir

#----- make absolute dir ---#
data=`readlink -f $data`
indir=`readlink -f $indir`

#----- convert PDB to XYZ files ---#
for i in `cat $data`
do
	$home/LIG_Tool/util/PDB_To_XYZ -i $indir/$i.pdb -a 1 -o $tmpdir/${i}_atom.xyz
done

#----- run pointsite --------#
source activate pointsite_inference
python $home/inference.py --output $outdir --data $tmpdir --select_list $data
source deactivate 2> /dev/null

#----- remove tmp -----#
rm -rf $tmpdir

#===== exit =====#
exit 0

