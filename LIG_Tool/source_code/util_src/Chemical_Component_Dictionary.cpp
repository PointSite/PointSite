#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <string.h>
using namespace std;


//----- remove 20 amino acid -----//
int Amino_Acid_Check(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 1;
	result=0;
	for(i=0;i<len;i++)
	{
		if(input[i]<'A' || input[i]>'Z')return 1;
		result+=(input[i]-'A')*(int)pow(26.0,1.0*i);
	}
	//switch
	switch(result)
	{
		case 286:return 0;
		case 4498:return 0;
		case 9256:return 0;
		case 10608:return 0;
		case 12794:return 0;
		case 9080:return 0;
		case 13812:return 0;
		case 16516:return 0;
		case 12383:return 0;
		case 2998:return 0;
		case 13635:return 0;
		case 12803:return 0;
		case 12960:return 0;
		case 2901:return 0;
		case 9921:return 0;
		case 11614:return 0;
		case 11693:return 0;
		case 10601:return 0;
		case 12135:return 0;
		case 7457:return 0;
		default:return 1;
	}
}




//====== given the Chemical Component Dictionary file, generate the corresponding mapping file =====//
//[note]: 
//[1] file is downloadable from "ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif"
//[2] in this file, only two fields are required, i) _chem_comp.id, and 2) _chem_comp.pdbx_type
//[3] according to CAMEO rule from "http://www.cameo3d.org/cameong_help/ligand-binding/"
//[3-1]  I(ions) = HETAI, HETIC, 
//[3-2]  O(organic) = HETAIN, ATOMS, HETAC, HETAD
//[3-3]  N(Poly-Nucleotides) = ATOMN 
//[3-4]  P(Peptides) = ATOMP 
//[4] the output format is as follows,
//[4-1] case 37951:return 0;
//[4-2] case 0:return "001";
//[4-3] case 0:return "HETAIN";
//[4-4] case 0:return "O";
//===============================// note over

//	temp=buf.substr(17,3); -> temp must have three letters !!!
void CCD_Head_Proc(string &temp,string &ws3)
{
	int ii;
	int len=temp.length();
	ws3.resize(3);
	for(ii=0;ii<3;ii++)ws3[ii]=' ';
	if(len>=3)
	{
		if(temp[0]==' ')//transfer
		{
			if(temp[1]==' ')
			{
				ws3[0]=temp[2];
			}
			else
			{
				ws3[0]=temp[1];
				ws3[1]=temp[2];
			}
		}
		else
		{
			for(ii=0;ii<3;ii++)ws3[ii]=temp[ii];
		}
	}
	else if(len==2)
	{
		if(temp[0]==' ')//transfer
		{
			ws3[0]=temp[1];
		}
		else
		{
			for(ii=0;ii<2;ii++)ws3[ii]=temp[ii];
		}
	}
	else if(len==1)
	{
		for(ii=0;ii<1;ii++)ws3[ii]=temp[ii];
	}
}
//---- ligand three digit transfer ---//
int Ligand_Trans(char code)
{
	switch(code)
	{
		case 'A': return 0;
		case 'B': return 1;
		case 'C': return 2;
		case 'D': return 3;
		case 'E': return 4;
		case 'F': return 5;
		case 'G': return 6;
		case 'H': return 7;
		case 'I': return 8;
		case 'J': return 9;
		case 'K': return 10;
		case 'L': return 11;
		case 'M': return 12;
		case 'N': return 13;
		case 'O': return 14;
		case 'P': return 15;
		case 'Q': return 16;
		case 'R': return 17;
		case 'S': return 18;
		case 'T': return 19;
		case 'U': return 20;
		case 'V': return 21;
		case 'W': return 22;
		case 'X': return 23;
		case 'Y': return 24;
		case 'Z': return 25;
		case '0': return 26;
		case '1': return 27;
		case '2': return 28;
		case '3': return 29;
		case '4': return 30;
		case '5': return 31;
		case '6': return 32;
		case '7': return 33;
		case '8': return 34;
		case '9': return 35;
		case ' ': return 36;
		case 'a': return 37;
		case 'b': return 38;
		case 'c': return 39;
		case 'd': return 40;
		case 'e': return 41;
		case 'f': return 42;
		case 'g': return 43;
		case 'h': return 44;
		case 'i': return 45;
		case 'j': return 46;
		case 'k': return 47;
		case 'l': return 48;
		case 'm': return 49;
		case 'n': return 50;
		case 'o': return 51;
		case 'p': return 52;
		case 'q': return 53;
		case 'r': return 54;
		case 's': return 55;
		case 't': return 56;
		case 'u': return 57;
		case 'v': return 58;
		case 'w': return 59;
		case 'x': return 60;
		case 'y': return 61;
		case 'z': return 62;
		default:return -1;
	}
}
//----- CAMEO four categories ----//
void CAMEO_Trans(string &input,string &output)
{
	if(input=="HETAI" || input=="HETIC")output="I";
	else if(input=="HETAIN" || input=="ATOMS" || input=="HETAC" || input=="HETAD")output="O";
	else if(input=="ATOMN")output="N";
	else if(input=="ATOMP")output="P";
	else output="X";
}

//------- main process ------//
int Chemical_Component_Dictionary_Parsing(string &infile,
	vector <string> &out1,vector <string> &out2,
	vector <string> &out3,vector <string> &out4)
{
	//read
	ifstream fin;
	string buf,temp;
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"infile %s not found!!\n",infile.c_str());
		exit(-1);
	}
	//load
	out1.clear();
	out2.clear();
	out3.clear();
	out4.clear();
	int count=0;
	int isHOH=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		stringstream www(buf);
		www>>temp;
		if(temp=="_chem_comp.id")
		{
			string value,retval;
			www>>value;
			CCD_Head_Proc(value,retval);
			if(retval=="HOH")isHOH=1;
			int retv=Amino_Acid_Check(retval.c_str());
			string outval=retval;
			if(retv==0)
			{
				//to lower case
				for(int i=0;i<(int)outval.length();i++)
					if(outval[i]>=65 && outval[i]<=90) outval[i]+=32;
			}
			//push to out2
			if(isHOH==0)
			{
				//get name
				stringstream sss1;
				sss1<<"case "<<count<<":return \""<<outval<<"\";";
				string finval1=sss1.str();
				out2.push_back(finval1);
				//get code
				int result=0;
				for(int i=0;i<3;i++)result+=(Ligand_Trans(outval[i]))*(int)pow(37.0,1.0*i);
				stringstream sss2;
				sss2<<"case "<<result<<":return "<<count<<";";
				string finval2=sss2.str();
				out1.push_back(finval2);
			}
		}
		if(temp=="_chem_comp.pdbx_type")
		{
			if(isHOH==0)
			{
				string value,retval;
				www>>value;
				CAMEO_Trans(value,retval);
				//get pdbx
				stringstream sss1;
				sss1<<"case "<<count<<":return \""<<value<<"\";";
				string finval1=sss1.str();
				out3.push_back(finval1);
				//get cameo
				stringstream sss2;
				sss2<<"case "<<count<<":return \""<<retval<<"\";";
				string finval2=sss2.str();
				out4.push_back(finval2);
				//count++
				count++;
			}
			isHOH=0;
		}
	}
	return count;
}

//---------- main ---------//
int main(int argc,char **argv)
{
	//---- Simp_Mapping ----//
	{
		if(argc<3)
		{
			fprintf(stderr,"Version 1.01 [Feb-19-2015] \n");
			fprintf(stderr,"Chemical_Component_Dictionary_Parse <ccd_file> <out_file> \n");
			fprintf(stderr,"[note]: \n");
			fprintf(stderr,"    [1] file is downloadable from ftp://ftp.wwpdb.org/pub/pdb/data/monomers/components.cif \n");
			fprintf(stderr,"    [2] in this file, only two fields are required, i) _chem_comp.id, and 2) _chem_comp.pdbx_type \n");
			fprintf(stderr,"    [3] according to CAMEO rule from http://www.cameo3d.org/cameong_help/ligand-binding/ \n");
			fprintf(stderr,"    [3-1]  I(ions) = HETAI, HETIC \n");
			fprintf(stderr,"    [3-2]  O(organic) = HETAIN, ATOMS, HETAC, HETAD \n");
			fprintf(stderr,"    [3-3]  N(Poly-Nucleotides) = ATOMN \n");
			fprintf(stderr,"    [3-4]  P(Peptides) = ATOMP \n");
			fprintf(stderr,"    [4] the output format is as follows, \n");
			fprintf(stderr,"    [4-1] case 37951:return 0; \n");
			fprintf(stderr,"    [4-2] case 0:return \"001\"; \n");
			fprintf(stderr,"    [4-3] case 0:return \"HETAIN\"; \n");
			exit(-1);
		}
		string ccd_file=argv[1];
		string out_file=argv[2];
		//process
		vector <string> out1;
		vector <string> out2;
		vector <string> out3;
		vector <string> out4;
		int count=Chemical_Component_Dictionary_Parsing(ccd_file,out1,out2,out3,out4);
		//output
		FILE *fp=fopen(out_file.c_str(),"wb");
		fprintf(fp,"#include \"Ligand_Utility.h\"\n");
		//out1
		fprintf(fp,"int Ligand_Num_Code(int code)\n");
		fprintf(fp,"{\n");
		fprintf(fp,"switch(code)\n");
		fprintf(fp,"{\n");
		for(int i=0;i<count;i++)fprintf(fp,"%s\n",out1[i].c_str());
		fprintf(fp,"default:return -1;\n");
		fprintf(fp,"}\n");
		fprintf(fp,"}\n");
		//out2
		fprintf(fp,"const char* Ligand_Num_Decode(int code)\n");
		fprintf(fp,"{\n");
		fprintf(fp,"switch(code)\n");
		fprintf(fp,"{\n");
		for(int i=0;i<count;i++)fprintf(fp,"%s\n",out2[i].c_str());
		fprintf(fp,"default:return \"UNK\";\n");
		fprintf(fp,"}\n");
		fprintf(fp,"}\n");
		//out3
		fprintf(fp,"const char* Ligand_PDBX(int code)\n");
		fprintf(fp,"{\n");
		fprintf(fp,"switch(code)\n");
		fprintf(fp,"{\n");
		for(int i=0;i<count;i++)fprintf(fp,"%s\n",out3[i].c_str());
		fprintf(fp,"default:return \"?\";\n");
		fprintf(fp,"}\n");
		fprintf(fp,"}\n");
		//out4
		fprintf(fp,"const char* Ligand_CAMEO(int code)\n");
		fprintf(fp,"{\n");
		fprintf(fp,"switch(code)\n");
		fprintf(fp,"{\n");
		for(int i=0;i<count;i++)fprintf(fp,"%s\n",out4[i].c_str());
		fprintf(fp,"default:return \"X\";\n");
		fprintf(fp,"}\n");
		fprintf(fp,"}\n");
		//fclose
		fclose(fp);
		//exit
		exit(0);
	}
}

