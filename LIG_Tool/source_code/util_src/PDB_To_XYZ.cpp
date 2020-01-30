#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <getopt.h>
using namespace std;


//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}

//------- int to string --------//
template <typename T>
string NumberToString( T Number )
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}

//------- remove space in a string -----//
void Kill_Space(string &in,string &out)
{
	int i;
	out="";
	for(i=0;i<(int)in.length();i++)
	{
		if(in[i]!=' ')out+=in[i];
	}
}
//----- check insertion code ------//
int Check_Ins(string &in)
{
	int i=(int)in.length()-1;
	if(in[i]>='0'&&in[i]<='9')return 0;
	else return 1;
}
//----- check chain_code ------//
int Check_Chain(string &in)
{
	int i=0;
	if(in[i]=='|')return 0;  //-> null chain
	else return 1;
}


//----------- SideChain_Related ------------//
//function (input three-digit atom)
int WWW_PDB_atom_name_hashing(const char *atom)
{
	int i;
	int result=0;
	int pos;
	int len=(int)strlen(atom);
	int rellen;
	// modify by wangsheng, 2011.04.30
	rellen=len<3?len:3;
	// modify by shaomingfu, 2011.04.12
	if(len!=3)len=3;
	for(i=0;i<rellen;i++)
	{
		if(atom[i]==' ')pos=0;
		else if(atom[i]>='A'&&atom[i]<='Z')pos=atom[i]-'A';
		else if(atom[i]>='0'&&atom[i]<='9')pos=atom[i]-'0';
		else pos=0;
		result+=(int)(pow(26.0,1.0*(len-i-1)))*pos;
	}
	return result;
}
//backbone 
int WWW_backbone_atom_name_encode(const char *atom)
{
	int pos=WWW_PDB_atom_name_hashing(atom);
	switch(pos)
	{
		case 8788:return 0;
		case 1352:
		{
			int len=(int)strlen(atom);
			if(len<2)return -1;
			char add=atom[1]; //second atom
			switch(add)
			{
				case 'A':return 1;
				case ' ':return 2;
				default:return -1;
			}
		}
		case 9464:return 3;
		default:return -1;
	}
}
//sidechain (ARNDCQEGHILKMFPSTWYVZ) 
int WWW_sidechain_atom_name_encode(const char *atom, char amino)
{
	int pos = WWW_PDB_atom_name_hashing(atom);
	switch(amino)
	{
		case 'A':
		{
			switch(pos)
			{
				case 1378:return 0;
				default:return -1;
			}
		}
		case 'R':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 8892:return 3;
				case 2002:return 4;
				case 8971:return 5;
				case 8972:return 6;
				default:return -1;
			}
		}
		case 'N':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 9543:return 2;
				case 8868:return 3;
				default:return -1;
			}
		}
		case 'D':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 9543:return 2;
				case 9544:return 3;
				default:return -1;
			}
		}
		case 'C':
		{
			switch(pos)
			{
				case 1378 :return 0;
				case 12324:return 1;
				default:return -1;
			}
		}
		case 'Q':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 9569:return 3;
				case 8894:return 4;
				default:return -1;
			}
		}
		case 'E':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 9569:return 3;
				case 9570:return 4;
				default:return -1;
			}
		}
		case 'G':
		{
			switch(pos)
			{
				case 1378:return 0; //note: pseudo CB(=CA)
				default:return -1;
			}
		}
		case 'H':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 8867:return 2;
				case 1432:return 3;
				case 1457:return 4;
				case 8894:return 5;
				default:return -1;
			}
		}
		case 'I':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1509:return 1;
				case 1510:return 2;
				case 1431:return 3;
				default:return -1;
			}
		}
		case 'L':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				default:return -1;
			}
		}
		case 'K':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				case 1456:return 3;
				case 9438:return 4;
				default:return -1;
			}
		}
		case 'M':
		{
			switch(pos)
			{
				case 1378 :return 0;
				case 1508 :return 1;
				case 12246:return 2;
				case 1456 :return 3;
				default:return -1;
			}
		}
		case 'F':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				case 1457:return 4;
				case 1458:return 5;
				case 2002:return 6;
				default:return -1;
			}
		}
		case 'P':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1430:return 2;
				default:return -1;
			}
		}
		case 'S':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 9620:return 1;
				default:return -1;
			}
		}
		case 'T':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 9621:return 1;
				case 1510:return 2;
				default:return -1;
			}
		}
		case 'W':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				case 8893:return 4;
				case 1458:return 5;
				case 1459:return 6;
				case 2004:return 7;
				case 2005:return 8;
				case 1536:return 9;
				default:return -1;
			}
		}
		case 'Y':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1508:return 1;
				case 1431:return 2;
				case 1432:return 3;
				case 1457:return 4;
				case 1458:return 5;
				case 2002:return 6;
				case 9646:return 7;
				default:return -1;
			}
		}
		case 'V':
		{
			switch(pos)
			{
				case 1378:return 0;
				case 1509:return 1;
				case 1510:return 2;
				default:return -1;
			}
		}
		case 'X':  //only consider CB!! (alanine model)
		{
			switch(pos)
			{
				case 1378:return 0;
				default:return -1;
			}
		}
		default:return -1;
	}
}


//----------- AA three to one ---------//
char WWW_Three2One_III(const char *input)
{
	int i;
	int len;
	int result;
	//encoding
	len=(int)strlen(input);
	if(len!=3)return 'X';
	result=0;
	for(i=0;i<len;i++)result+=(input[i]-'A')*(int)pow(1.0*26,1.0*i);
	//switch
	switch(result)
	{
		case 286:return 'A';
		case 4498:return 'R';
		case 9256:return 'N';
		case 10608:return 'D';
		case 12794:return 'C';
		case 9080:return 'Q';
		case 13812:return 'E';
		case 16516:return 'G';
		case 12383:return 'H';
		case 2998:return 'I';
		case 13635:return 'L';
		case 12803:return 'K';
		case 12960:return 'M';
		case 2901:return 'F';
		case 9921:return 'P';
		case 11614:return 'S';
		case 11693:return 'T';
		case 10601:return 'W';
		case 12135:return 'Y';
		case 7457:return 'V';
		default:return 'X';
	}
}

//========= XYZ format for point-cloud ==========//
/*
    6 DOg  -22.389   29.406   56.176
    6 DOh  -23.720   31.160   56.319
    7 KNa  -22.679   27.076   60.272
    7 KCb  -21.836   26.578   61.340
    7 KCc  -20.356   26.635   60.950
    7 KOd  -19.746   25.631   60.582
    7 KCe  -22.307   25.189   61.783
    7 KCf  -23.751   25.237   62.298
    7 KCg  -24.043   24.207   63.371
    7 KCh  -25.357   24.538   64.095
    7 KNi  -25.610   23.600   65.237
    8 QNa  -19.785   27.829   61.073
    8 QCb  -18.384   28.110   60.732
    8 QCc  -17.293   27.383   61.502
    8 QOd  -16.145   27.388   61.063
...
*/

//--------- PDB_To_XYZ ----------//
//-> PDB format:
/*
ATOM    496  N   MET A   1     -30.766  33.288  73.186  1.00 19.31      A    N
ATOM    497  CA  MET A   1     -29.721  33.119  72.137  1.00 18.15      A    C
ATOM    498  C   MET A   1     -30.322  33.150  70.741  1.00 19.44      A    C
ATOM    499  O   MET A   1     -31.535  33.139  70.580  1.00 20.69      A    O
ATOM    500  CB  MET A   1     -28.881  31.843  72.350  1.00 16.43      A    C
ATOM    501  CG  MET A   1     -29.636  30.528  72.448  1.00 17.15      A    C
...
*/
void PDB_To_XYZ(string &pdb,FILE *fp,
	int resi_start,int resi_fix,int addi_columns)
{
	ifstream fin;
	string buf,temp,name;
	//read
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",pdb.c_str());
		exit(-1);
	}
	int len;
	int resi;
	string type;
	string prev="";
	int first=1;
	int resi_num=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<3)continue;
		//check TER
		temp=buf.substr(0,3);
		if(temp=="TER"||temp=="END")break;
		//check ATOM
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp!="ATOM"&&temp!="HETA")continue;
		type=temp;
		//ami
		temp=buf.substr(17,3);
		char c=WWW_Three2One_III(temp.c_str());
		char h=buf[77];
		//chain
		char chain=buf[21];
		//resi
		string curr=buf.substr(22,5);
		if(first==1)
		{
			first=0;
			prev=curr;
		}
		if(curr!=prev)
		{
			resi_num++;
			prev=curr;
		}
		string resi_str;
		if(resi_start<0 && resi_fix<0)
		{
			resi_str=curr;
		}
		else
		{
			if(resi_fix>=0)  //-> fix residue (has priority over resi_start)
			{
				resi_str=NumberToString(resi_fix);
			}
			else             //-> use resi_start
			{
				resi_str=NumberToString(resi_start+resi_num);
			}
		}
		//-> reformat to A|123 style //-> start
		string resi_str_="";
		resi_str_.push_back(chain);
		resi_str_=resi_str_+"|"+resi_str;
		Kill_Space(resi_str_,resi_str);
		if(Check_Ins(resi_str)!=1)resi_str.push_back(' ');
		if(Check_Chain(resi_str)==0)resi_str=" "+resi_str;
		//atom
		char atom_char='!';
		if(type=="ATOM"&&c!='X')
		{
			temp=buf.substr(13,3);
			int atom_encode=0;
			int back_ret=WWW_backbone_atom_name_encode(temp.c_str());
			int side_ret=WWW_sidechain_atom_name_encode(temp.c_str(),c);
			if((back_ret<0&&side_ret<0)||(back_ret>=0&&side_ret>=0)) //BAD!!
			{
				fprintf(stderr,"WARNING!!! BAD_ATOM!!!\n%s\n",buf.c_str());
			}
			else
			{
				if(back_ret>=0)atom_encode=back_ret;
				if(side_ret>=0)atom_encode=4+side_ret;
			}
			atom_char=atom_encode+'a';
		}
		//mol
		temp=buf.substr(30,8);
		double x=atof(temp.c_str());
		temp=buf.substr(38,8);
		double y=atof(temp.c_str());
		temp=buf.substr(46,8);
		double z=atof(temp.c_str());
		//bfactor
		temp=buf.substr(60,6);
		double bfactor=atof(temp.c_str());
		//output
		if(addi_columns<0)
		{
			fprintf(fp,"%7s %c%c%c %8.3f %8.3f %8.3f \n",
				resi_str.c_str(),c,h,atom_char,x,y,z);
		}
		else
		{
			fprintf(fp,"%7s %c%c%c %8.3f %8.3f %8.3f 0 %5.2f \n",
				resi_str.c_str(),c,h,atom_char,x,y,z,bfactor);
		}
	}
}


//---------- usage ---------//
void Usage()
{
	fprintf(stderr,"Version: 0.98 \n");
	fprintf(stderr,"PDB_To_XYZ -i pdb_input -o xyz_output [-s resi_start] [-f resi_fix] [-a addi_col] \n\n");
	fprintf(stderr,"Usage : \n\n");
	fprintf(stderr,"-i pdb_input :         Input PDB file. \n\n");
	fprintf(stderr,"-o xyz_output :        Output XYZ file. \n\n");
	fprintf(stderr,"-s resi_start :        Default: set -1 to use original residue numbering, \n");
	fprintf(stderr,"                       or, set the starting residue (1-base). \n\n");
	fprintf(stderr,"-f resi_fix :          Default: set -1 to use original residue numbering, \n");
	fprintf(stderr,"                       or, fix the residue number to the given value. \n\n");
	fprintf(stderr,"-a addi_col :          Default: set -1 to CANCEL output additional columns, \n");
	fprintf(stderr,"                       or, specify 1 to output '0 bfactor ' columns. \n\n");
}


//-------- main ---------//
int main(int argc,char **argv)
{
	//------ PDB_To_XYZ -------//
	{
		if(argc<3)
		{
			Usage();
			exit(-1);
		}
		string pdb_input="";
		string xyz_output="";
		int resi_start=-1;
		int resi_fix=-1;
		int addi_col=-1;
		//command-line arguments process
		extern char* optarg;
		char c = 0;
		while ((c = getopt(argc, argv, "i:o:s:f:a:")) != EOF) 
		{
			switch (c) 
			{
				case 'i':
					pdb_input = optarg;
					break;
				case 'o':
					xyz_output = optarg;
					break;
				case 's':
					resi_start = atoi(optarg);
					break;
				case 'f':
					resi_fix = atoi(optarg);
					break;
				case 'a':
					addi_col = atoi(optarg);
					break;
				default:
					Usage();
					exit(-1);
			}
		}
		//check arguments
		if(pdb_input=="")
		{
			fprintf(stderr,"pdb_input should be specified \n");
			exit(-1);
		}
		if(xyz_output=="")
		{
			fprintf(stderr,"xyz_output should be specified \n");
			exit(-1);
		}
		//proc
		FILE *fp=fopen(xyz_output.c_str(),"wb");
		PDB_To_XYZ(pdb_input,fp,resi_start,resi_fix,addi_col);
		fclose(fp);
		//exit
		exit(0);
	}
}



