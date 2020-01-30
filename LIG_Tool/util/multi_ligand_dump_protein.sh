#!/bin/bash

# ----- usage ------ #
function usage()
{
	echo "multi_ligand_dump (protein-level) v0.10 [Sep-07-2019] "
	echo "    A simple script to dump '0,1,2,..' ligand labels in XYZ format"
	echo ""
	echo "USAGE:  ./multi_ligand_dump.sh <-i in_list> <-o out_root> <-r lig_root> " 
	echo "                   [-d distance_cut] [-H home] "
	echo "Options:"
	echo ""
	echo "***** required arguments *****"
	echo "-i in_list        : Input PDB list in 1pdb format. "
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
	#--- check XYZ files in the lig/ folder ---#
	has_xyz=1;
	for k in `awk '{print $1}' $lig/$i.chain_log`;
	do
		if [ ! -s "$lig/$k.pc_xyz" ]  #-> use already exist pc_xyz file
		then
			has_xyz=0;
		fi
		cat $lig/$k.pc_xyz >> $i.xyz_protein;
	done;

	#--- cat all chains into a complex ---#
	if [ $has_xyz -eq 1 ]
	then
		rm -f $i.xyz_protein;
		for k in `awk '{print $1}' $lig/$i.chain_log`;
		do
			cat $lig/$k.pc_xyz >> $i.xyz_protein;
		done;
	else
		rm -f $i.pdb_protein;
		for k in `awk '{print $1}' $lig/$i.chain_log`;
		do
			cat $lig/$k.pdb >> $i.pdb_protein;
		done;
		$home/PDB_To_XYZ -i $i.pdb_protein -o $i.xyz_protein -a 1;
		rm -f $i.pdb_protein;
	fi

	#--- cat all ligands into a complex ---#
	rm -f $i.xyz_lig;
	rm -f $i.ligand;
	count=1;
	for k in `awk '{print $1}' $lig/$i.ligand_size`;
	do
		b=${k:5:9};
		$home/PDB_To_XYZ -i $lig/$k.pdb -o $k.xyz -f $count;
		((count++));
		cat $k.xyz >> $i.xyz_lig;
		echo $b >> $i.ligand;
		rm -f $k.xyz;
	done;

	#--- call XYZ_ContResi -----#
	$home/XYZ_ContResi $i.xyz_protein $i.xyz_lig $distance_cut $out/${i}_atom.xyz 1 > $out/${i}_resi;
	rm -f $i.xyz_protein;
	mv $i.xyz_lig $out/${i}_lig.xyz;
	mv $i.ligand $out;
done


