CXX = g++
#CCFLAGS = -O3 -Wall -ffast-math -static
CCFLAGS = -O3
OBJ = Utility.o Computation_Utility.o Confo_Back.o Confo_Beta.o Confo_Lett.o Acc_Surface.o  getopt.o Hydro_Bond.o LIG_Tool.o Mol_File.o Mol_Ori.o Mol_Out.o PDB_Utility.o Backbone_Sidechain.o PDB_Residue.o PDB_Chain.o PDB_Chain_Fold.o PDB_File.o PhiPsi_Trans.o XYZ_DBT.o XYZ.o Ligand_Utility.o

make: $(OBJ)
	${CXX} $(CCFLAGS)  -o LIG_Tool $(OBJ)
#	rm -f $(OBJ)

Utility.o  : Utility.h
	${CXX} ${CCFLAGS} -c Utility.cpp
Computation_Utility.o  : Computation_Utility.h
	${CXX} ${CCFLAGS} -c Computation_Utility.cpp
Confo_Back.o  : Confo_Back.h
	${CXX} ${CCFLAGS} -c Confo_Back.cpp
Confo_Beta.o  : Confo_Beta.h
	${CXX} ${CCFLAGS} -c Confo_Beta.cpp
Confo_Lett.o  : Confo_Lett.h
	${CXX} ${CCFLAGS} -c Confo_Lett.cpp
Acc_Surface.o : Acc_Surface.h
	${CXX} ${CCFLAGS} -c Acc_Surface.cpp
getopt.o  : getopt.h
	${CXX} ${CCFLAGS} -c getopt.cpp
Hydro_Bond.o  : Hydro_Bond.h
	${CXX} ${CCFLAGS} -c Hydro_Bond.cpp
LIG_Tool.o  : 
	${CXX} ${CCFLAGS} -c LIG_Tool.cpp
Mol_File.o  : Mol_File.h
	${CXX} ${CCFLAGS} -c Mol_File.cpp
Mol_Ori.o  : Mol_Ori.h
	${CXX} ${CCFLAGS} -c Mol_Ori.cpp
Mol_Out.o  : Mol_Out.h
	${CXX} ${CCFLAGS} -c Mol_Out.cpp	
PDB_Utility.o  : PDB_Utility.h
	${CXX} ${CCFLAGS} -c PDB_Utility.cpp
Backbone_Sidechain.o  : Backbone_Sidechain.h
	${CXX} ${CCFLAGS} -c Backbone_Sidechain.cpp
PDB_Residue.o  : PDB_Residue.h
	${CXX} ${CCFLAGS} -c PDB_Residue.cpp
PDB_Chain.o  : PDB_Chain.h
	${CXX} ${CCFLAGS} -c PDB_Chain.cpp
PDB_Chain_Fold.o  : PDB_Chain_Fold.h
	${CXX} ${CCFLAGS} -c PDB_Chain_Fold.cpp	
PDB_File.o  : PDB_File.h
	${CXX} ${CCFLAGS} -c PDB_File.cpp
PhiPsi_Trans.o  : PhiPsi_Trans.h
	${CXX} ${CCFLAGS} -c PhiPsi_Trans.cpp
XYZ_DBT.o  : XYZ_DBT.h
	${CXX} ${CCFLAGS} -c XYZ_DBT.cpp
XYZ.o  : XYZ.h
	${CXX} ${CCFLAGS} -c XYZ.cpp
Ligand_Utility.o : Ligand_Utility.h
	${CXX} ${CCFLAGS} -c Ligand_Utility.cpp

clean:
	rm -f LIG_Tool $(OBJ)
