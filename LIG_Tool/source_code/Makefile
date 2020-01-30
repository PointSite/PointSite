

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

LIGToolDir = ./LIG_Tool_src
LIGToolEXE = /LIG_Tool
UTILS = ./util_src/

all: LIG_Tool util

LIG_Tool:
	cd $(LIGToolDir)&&make
	mv ${LIGToolDir}$(LIGToolEXE) ../
	rm -f ${LIGToolDir}/*.o

util:
	$(CXX) $(UTILS)Chemical_Component_Dictionary.cpp $(FLAG) -o Chemical_Component_Dictionary
	$(CXX) $(UTILS)XYZ_ContResi.cpp $(FLAG) -o XYZ_ContResi
	$(CXX) $(UTILS)PDB_To_XYZ.cpp $(FLAG) -o PDB_To_XYZ
	$(CXX) $(UTILS)XYZ_To_PDB.cpp $(FLAG) -o XYZ_To_PDB
	$(CXX) $(UTILS)XYZ_Gen_Feature.cpp $(FLAG) -o XYZ_Gen_Feature
	$(CXX) $(UTILS)XYZ_DCCnDCA.cpp $(FLAG) -o XYZ_DCCnDCA
	mkdir -p ../util
	mv Chemical_Component_Dictionary XYZ_ContResi PDB_To_XYZ XYZ_To_PDB XYZ_Gen_Feature XYZ_DCCnDCA ../util/

# clean all the temp files
clean:
