#!/bin/bash

# ----- usage ------ #
function usage()
{
	echo "multi_ligand_dump (chain-level) v0.10 [Sep-07-2019] "
	echo "    A simple script to dump '0,1,2,..' ligand labels in XYZ format"
	echo ""
	echo "USAGE:  ./multi_ligand_dump.sh <-i in_list> <-o out_root> <-r lig_root> " 
	echo "                   [-d distance_cut] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i in_list        : Input PDB list in 1pdbA format. "
	echo ""
	echo "-o out_root       : Output root for dumped XYZ files. "
	echo ""
	echo "-r lig_root       : The root containing processed ligands. "
	echo ""
	echo "***** optional arguments (database) *****"
	echo "-d distance_cut   : ligand distance cutoff [default: 6.5]  "
	echo ""
	echo "***** home directory *****"
	echo "-H home           : home directory of multi_ligand_dump.sh "
	echo "                    [default = `dirname $0`] "
	echo ""
	exit 1
}


#------------------------------------------------------------#
##### ===== get pwd and check BlastSearchHome ====== #########
#------------------------------------------------------------#

#------ current directory ------#
curdir="$(pwd)"

#-------- check usage -------#
if [ $# -lt 1 ];
then
	usage
fi


#---------------------------------------------------------#
##### ===== All arguments are defined here ====== #########
#---------------------------------------------------------#

#------- required arguments ------------#
#-> input output
in_list=""          #-> input PDB list in 1pdbA format
out_root=""         #-> output root for dumped XYZ files
lig_root=""         #-> ligand root after executing LIG_Tool
#-> parameters
distance_cut=6.5    #-> 6.5 is the cutoff used in scPDB
#-> home
home=`dirname $0`   #-> home directory

#------- parse arguments ---------------#
while getopts ":i:o:r:d:H:" opt;
do
	case $opt in
	#-> required arguments
	i)
		in_list=$OPTARG
		;;
	o)
		out_root=$OPTARG
		;;
	r)
		lig_root=$OPTARG
		;;
	#-> parameters
	d)
		distance_cut=$OPTARG
		;;
	#-> home directory
	H)
		home=$OPTARG
		;;
	#-> default
	\?)
		echo "Invalid option: -$OPTARG" >&2
		exit 1
		;;
	:)
		echo "Option -$OPTARG requires an argument." >&2
		exit 1
		;;
	esac
done


#---------------------------------------------------------#
##### ===== Part 0: initial argument check ====== #########
#---------------------------------------------------------#

# ------ check home directory ---------- #
if [ ! -d "$home" ]
then
	echo "home directory $home not exist " >&2
	exit 1
fi
home=`readlink -f $home`

# ------ check input list ------#
if [ ! -s "$in_list" ]
then
	echo "in_list $in_list not found !!" >&2
	exit 1
fi
in_list=`readlink -f $in_list`
fulnam=`basename $in_list`
relnam=${fulnam%.*}

# ------ check output root ------#
if [ "$out_root" == "" ]
then
	out_root=${relnam}.multi_ligand_dump
fi
mkdir -p $out_root
out_root=`readlink -f $out_root`

# ------ check ligand root ------#
if [ ! -d "$lig_root" ]
then
	echo "lig_root $lig_root not found !!" >&2
	exit 1
fi
lig_root=`readlink -f $lig_root`


#-----------------------------------------------------------#
##### ===== Part 1: run multi_ligand_dump.sh ====== #########
#-----------------------------------------------------------#

#---- preliminary ---#
#-> lig root
lig=$lig_root               #e.g.,   ../lig
#-> input/output
list=$in_list               #e.g.,   multilig_less3_list
out=$out_root               #e.g.,   multilig_less3_xyz
mkdir -p $out

#---- process -------#
for i in `cat $list`;
do
	id=${i:0:4};
	rm -f $i.xyz_lig;
	rm -f $i.ligand;
	count=1;
	for k in `grep \^${i} $lig/$id.ligand_log | awk '{print $1}'`;
	do
		a=${k:0:4};
		b=${k:5:9};
		$home/PDB_To_XYZ -i $lig/"$a$b".pdb -o $i.xyz_ -f $count;
		((count++));
		cat $i.xyz_ >> $i.xyz_lig;
		echo $b >> $i.ligand;
		rm -f $i.xyz_;
	done;
	if [ -s "$lig/$i.pc_xyz" ]  #-> use already exist pc_xyz file
	then
		$home/XYZ_ContResi $lig/$i.pc_xyz $i.xyz_lig $distance_cut $out/${i}_atom.xyz 1 > $out/${i}_resi;
	else                        #-> create xyz file from PDB file
		$home/PDB_To_XYZ -i $lig/$i.pdb -o $i.xyz_atom -a 1;
		$home/XYZ_ContResi $i.xyz_atom $i.xyz_lig $distance_cut $out/${i}_atom.xyz 1 > $out/${i}_resi;
	fi
	rm -f $i.xyz_atom;
	mv $i.xyz_lig $out/${i}_lig.xyz;
	mv $i.ligand $out;
done


