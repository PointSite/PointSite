#!/bin/bash

if [ $# -lt 1 ]
then
	echo "Usage: ./update_Ligand_Utility.sh <ANYKEY> "
	echo "[note]: if ANYKEY is set to -1, then remove temporary folder which contains the file 'components.cif'"
	exit 1
fi

#part0 -> create a temporary folder
DATE=`date '+%Y_%m_%d'`
tmp="update_Ligand_Utility_${DATE}"
mkdir -p $tmp

#part1 -> download components.cif from ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif
rm -f $tmp/components.cif
wget -q ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif.gz -O $tmp/components.cif.gz
gunzip -q $tmp/components.cif.gz 
if [ ! -f "$tmp/components.cif" ]
then
	echo "file components.cif download error !! "
	exit 1
fi
echo "part1 done -> download components.cif"

#part2 -> generate Ligand_Utility.cpp from components.cif by ./Chemical_Component_Dictionary
cp source_code/Ligand_Utility_cpp/Ligand_Utility.cpp source_code/Ligand_Utility_cpp/Ligand_Utility.cpp_$DATE
util/Chemical_Component_Dictionary $tmp/components.cif $tmp/Ligand_Utility.cpp
cp $tmp/Ligand_Utility.cpp source_code/Ligand_Utility_cpp/
echo "part2 done -> generate Ligand_Utility.cpp"

#part3 -> remove temporary folder
if [ "$1" == "-1" ]
then
	rm -rf $tmp
fi

#-------- exit --------#
exit 0

