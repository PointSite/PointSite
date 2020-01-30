#include <map>
#include "PDB_Chain_Fold.h"
#include "Confo_Back.h"
#include "Confo_Beta.h"
#include "Confo_Lett.h"
#include "Acc_Surface.h"
#include "Mol_File.h"
#include "Mol_Out.h"
#include "getopt.h"
#include "Ligand_Utility.h"
#include "PDB_Utility.h"
using namespace std;


//-----------------------------------------------------------------------------------------------------------//
//---- print_help_msg => print help message
void print_help_msg(void) 
{
	cout << "========================================================|" << endl;
	cout << "LIG_Tool  (version 1.06) [2019.10.03]                   |" << endl;
	cout << "    Extract ligands and chains from official PDB file   |" << endl;
	cout << "Usage: ./LIG_Tool <-i input_pdb> [-o out_name]          |" << endl;
	cout << "       [-p chain_root] [-q ligand_root] [-L log_root]   |" << endl;
	cout << "       [-n length_cut] [-d distance_cut] [-l log]       |" << endl;
	cout << "       [-O PC_out] [-T PC_type] [-t LG_type]            |" << endl;
	cout << "       [-f filter] [-m mincut] [-M atomcut] [-N minnum] |" << endl;
	cout << "--------------------------------------------------------|" << endl;
	cout << "-i input_pdb : input original PDB file, e.g., 1col.pdb  |" << endl;
	cout << "-o out_name  : output file name. [default: input_name]  |" << endl;
	cout << "-p chain_root  : chain output root [default: ./ ]       |" << endl;
	cout << "-q ligand_root : ligand output root [default: ./ ]      |" << endl;
	cout << "-n length_cut  : chain length cutoff [default: 25 ]     |" << endl;
	cout << "-d distance_cut: ligand distance cutoff [default: 4.0]  |" << endl;
	cout << "-l log : output log files (1) or not (0) [default: 0]   |" << endl;
	cout << "-L log_root: log output root [default:null]             |" << endl;
	cout << "-O PC_root : binding residues in PC [default:null]      |" << endl;
	cout << "-T PC_type : CA+CB (-1), backbone+CB [0], full atom (1) |" << endl;
	cout << "-t LG_type : select ligand type. [ default: 'IOPNX']    |" << endl;
	cout << "-f filter  : list for filtered ligands. [default:null]  |" << endl;
	cout << "-m mincut  : min_num of binding residues. [default:1]   |" << endl;
	cout << "-M atomcut : min_num of binding atoms. [default:1]      |" << endl;
	cout << "-N minnum  : min_num of the ligand atoms. [default:1]   |" << endl;
	cout << "========================================================|" << endl;
	exit(-1);
}
//---- WebServer's default input ----//
string INPUT_FILE="";
string OUTPUT_NAME="";
string CHAIN_OUTROOT="./";  //chain output root
string LIGAND_OUTROOT="./"; //ligand output root
int LENGTH_CUT=25;          //chain length cutoff
double DISTANCE_CUT=4.0;    //ligand distance cutoff
int LOG_OR_NOT=0;           //default: don't output log
string LOG_ROOT="";         //log output root
string PC_ROOT="";          //point cloud output root
int PC_TYPE=0;              //default: backbone+CB
string LG_TYPE="IOPNX";     //default: ALL possible ligand type
string LG_FILTER="";        //filtered ligands (default: null)
int LG_MINCUT=1;            //minimal number of ligand binding residues (default:1)
int LG_ATOMCUT=1;           //minimal number of ligand binding atoms (default:1)
int LG_MINNUM=1;            //minimal number of the ligand atoms (default:1)

//-----------------------------------------------------------------------------------------------------------//
//---- parameter editor ----//
static option long_options[] =
{
	{"input",   required_argument, NULL, 'i'},
	{"output",  no_argument,       NULL, 'o'},
	{"chain",   no_argument,       NULL, 'p'},
	{"ligand",  no_argument,       NULL, 'q'},
	{"length",  no_argument,       NULL, 'n'},
	{"dist",    no_argument,       NULL, 'd'},
	{"log",     no_argument,       NULL, 'l'},
	{"LOGroot", no_argument,       NULL, 'L'},
	{"PCroot",  no_argument,       NULL, 'O'},
	{"PCtype",  no_argument,       NULL, 'T'},
	{"LGtype",  no_argument,       NULL, 't'},
	{"LGfilt",  no_argument,       NULL, 'f'},
	{"LGmin",   no_argument,       NULL, 'm'},
	{"ATmin",   no_argument,       NULL, 'M'},
	{"LGnum",   no_argument,       NULL, 'N'},
	{0, 0, 0, 0}
};
//-----------------------------------------------------------------------------------------------------------//
//---- process_args => process input parameter args
void process_args(int argc,char** argv) 
{
	string buf;
	int opt;
	if(1==argc)print_help_msg();    
	while(true) 
	{
		int option_index=0;
		opt=getopt_long(argc,argv,"i:o:p:q:n:d:l:L:O:T:t:f:m:M:N:",
			   long_options,&option_index);
		if (opt==-1)break;	
		switch(opt) 
		{
			case 'i':
				INPUT_FILE=optarg;
				break;
			case 'o':
				OUTPUT_NAME=optarg;
				break;
			case 'p':
				CHAIN_OUTROOT=optarg;
				break;
			case 'q':
				LIGAND_OUTROOT=optarg;
				break;
			case 'n':
				LENGTH_CUT=atoi(optarg);
				break;
			case 'd':
				DISTANCE_CUT=atof(optarg);
				break;
			case 'l':
				LOG_OR_NOT=atoi(optarg);
				break;
			case 'L':
				LOG_ROOT=optarg;
				break;
			case 'O':
				PC_ROOT=optarg;
				break;
			case 'T':
				PC_TYPE=atoi(optarg);
				break;
			case 't':
				LG_TYPE=optarg;
				break;
			case 'f':
				LG_FILTER=optarg;
				break;
			case 'm':
				LG_MINCUT=atoi(optarg);
				break;
			case 'M':
				LG_ATOMCUT=atoi(optarg);
				break;
			case 'N':
				LG_MINNUM=atoi(optarg);
				break;
			default:
				exit(-1);
		}
	}
}

//====== data structure ======//
struct Ligand_Struc
{
	string lig_name;
	string lig_type;
	int lig_moln;
	vector <XYZ> lig_xyz;
	vector <string> lig_data;
};

//--------- Ligand_Mapping --------//
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
		default:return -1; 
	} 
}

//--------- Ligand Type Mapping --------//
int Ligand_Type(char code) 
{ 
	switch(code) 
	{
		case 'I': return 0;
		case 'O': return 1;
		case 'P': return 2;
		case 'N': return 3;
		case 'X': return 4;
		default:return -1; 
	} 
}


//================== All_Ligands_Process =============//__110530__//
//-> ligand_string_to_xyz
int Ligand_String_To_XYZ(vector <string> &input, vector <XYZ> &output)
{
	int i;
	int size=(int)input.size();
	output.clear();
	string buf,temp;
	for(i=0;i<size;i++)
	{
		buf=input[i];
		temp=buf.substr(30,8);
		XYZ xyz;
		xyz.X=atof(temp.c_str());
		temp=buf.substr(38,8);
		xyz.Y=atof(temp.c_str());
		temp=buf.substr(46,8);
		xyz.Z=atof(temp.c_str());
		output.push_back(xyz);
	}
	return size;
}
//--------- extract ligand from PDB ------------//
int PDB_Extract_Ligand(string &pdb,vector <Ligand_Struc> &output)
{
	//---- get MODRES ----//
	Mol_File mol_input;
	map <string,string > modres_map;
	{
		const string sss=pdb;
		mol_input.Process_MODRES_Mapping(sss,modres_map);
	}
	//---- get MODRES ----//over

	output.clear();
	//--- list for mapping ---//
	map<string, int > in_mapping;
	map<string, int>::iterator iter;
	in_mapping.clear();
	//--- data for mapping ---//
	vector <int> in_mapp_int;
	vector <string> in_mapp_nam;
	vector <vector < vector <string> > > in_mapp_string;
	in_mapp_int.clear();
	in_mapp_nam.clear();
	in_mapp_string.clear();
	string lig_dummy;

	//--- list for mapping ---//
	ifstream fin;
	string buf,temp;
	string in_nam;
	getBaseName(pdb,in_nam,'/','.');
	//read
	fin.open(pdb.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",pdb.c_str());
		return -1;
	}
	int i;
	int ii;
	int len;
	int first=1;
	int tot_lig=0;
	int key;
	char ws3[4];
	char chain;
	vector <string> in_record;
	in_record.clear();
	int prev=-1;
	string prev_chain;
	string current_chain;
	char in_prev[4];
	//processs
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		len=(int)buf.length();
		if(len<4)continue;
		temp=buf.substr(0,4);
		if(temp=="END " || temp=="ENDM")break; //this might be modified in the future
		if(temp=="TER ") //output "ter"
		{
			//output
			if(first==0)
			{
				first=1;
				//record
				vector <string> in_record_tmp;
				for(i=0;i<(int)in_record.size();i++)
				{
					if(in_record[i][77]!='H')
					{
						in_record_tmp.push_back(in_record[i]);
					}
				}
				//output
				lig_dummy="";
				lig_dummy=lig_dummy+in_prev;
				iter = in_mapping.find(lig_dummy);
				if(iter != in_mapping.end()) //already exist ligand
				{
					key=in_mapping[lig_dummy];
					in_mapp_int[key-1]++;
					in_mapp_string[key-1].push_back(in_record_tmp);
				}
				else                         //newly appear ligand
				{
					tot_lig++;
					in_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
					string wsss=lig_dummy;
					for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
					in_mapp_nam.push_back(wsss);
					in_mapp_int.push_back(1);
					vector < vector <string> > dummy_vec;
					in_mapp_string.push_back(dummy_vec);
					in_mapp_string[tot_lig-1].push_back(in_record_tmp);
				}
			}
		}
		//get atom
		if(temp!="ATOM" && temp!="HETA")continue;
		if(len<20)continue;
		temp=buf.substr(17,3);
		if(temp[0]==' ')//transfer
		{
			for(ii=0;ii<3;ii++)ws3[ii]=' ';
			ws3[ii]='\0';
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
		else strcpy(ws3,temp.c_str());

		//MODRES MAP
		{
			string temp_ori=ws3;
			string temp_mod=temp_ori;
			mol_input.MODRES_Map(temp_ori,temp_mod,modres_map);
			for(ii=0;ii<3;ii++)ws3[ii]=temp_mod[ii];
		}
		//MODRES MAP over

		//record
		chain=buf[21];
		current_chain=buf.substr(21,6);
		int result=0;
		for(i=0;i<3;i++)result+=(Ligand_Trans(ws3[i]))*(int)pow(37.0,1.0*i);
		int retv=Ligand_Num_Code(result);
		if(retv!=-1)
		{
			if(first==1)
			{
				first=0;
				prev=retv;
				prev_chain=current_chain;
				strcpy(in_prev,ws3);
				in_record.clear();
				in_record.push_back(buf);
			}
			else
			{
				if(retv==prev && current_chain==prev_chain)in_record.push_back(buf);
				else
				{
					//record
					vector <string> in_record_tmp;
					for(i=0;i<(int)in_record.size();i++)
					{
						if(in_record[i][77]!='H')
						{
							in_record_tmp.push_back(in_record[i]);
						}
					}
					//output
					lig_dummy="";
					lig_dummy=lig_dummy+in_prev;
					iter = in_mapping.find(lig_dummy);
					if(iter != in_mapping.end()) //already exist ligand
					{
						key=in_mapping[lig_dummy];
						in_mapp_int[key-1]++;
						in_mapp_string[key-1].push_back(in_record_tmp);
					}
					else                         //newly appear ligand
					{
						tot_lig++;
						in_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
						string wsss=lig_dummy;
						for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
						in_mapp_nam.push_back(wsss);
						in_mapp_int.push_back(1);
						vector < vector <string> > dummy_vec;
						in_mapp_string.push_back(dummy_vec);
						in_mapp_string[tot_lig-1].push_back(in_record_tmp);
					}
					//continue
					prev=retv;
					prev_chain=current_chain;
					strcpy(in_prev,ws3);
					in_record.clear();
					in_record.push_back(buf);
				}
			}
		}
		else
		{
			if(first==0)
			{
				first=1;
				//record
				vector <string> in_record_tmp;
				for(i=0;i<(int)in_record.size();i++)
				{
					if(in_record[i][77]!='H')
					{
						in_record_tmp.push_back(in_record[i]);
					}
				}
				//output
				lig_dummy="";
				lig_dummy=lig_dummy+in_prev;
				iter = in_mapping.find(lig_dummy);
				if(iter != in_mapping.end()) //already exist ligand
				{
					key=in_mapping[lig_dummy];
					in_mapp_int[key-1]++;
					in_mapp_string[key-1].push_back(in_record_tmp);
				}
				else                         //newly appear ligand
				{
					tot_lig++;
					in_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
					string wsss=lig_dummy;
					for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
					in_mapp_nam.push_back(wsss);
					in_mapp_int.push_back(1);
					vector < vector <string> > dummy_vec;
					in_mapp_string.push_back(dummy_vec);
					in_mapp_string[tot_lig-1].push_back(in_record_tmp);
				}
			}
		}
	}
	//final
	if(first==0)
	{
		first=1;
		//record
		vector <string> in_record_tmp;
		for(i=0;i<(int)in_record.size();i++)
		{
			if(in_record[i][77]!='H')
			{
				in_record_tmp.push_back(in_record[i]);
			}
		}
		//output
		lig_dummy="";
		lig_dummy=lig_dummy+in_prev;
		iter = in_mapping.find(lig_dummy);
		if(iter != in_mapping.end()) //already exist ligand
		{
			key=in_mapping[lig_dummy];
			in_mapp_int[key-1]++;
			in_mapp_string[key-1].push_back(in_record_tmp);
		}
		else                         //newly appear ligand
		{
			tot_lig++;
			in_mapping.insert(map < string, int >::value_type(lig_dummy, tot_lig));
			string wsss=lig_dummy;
			for(ii=0;ii<(int)wsss.length();ii++)if(wsss[ii]==' ')wsss[ii]='_';
			in_mapp_nam.push_back(wsss);
			in_mapp_int.push_back(1);
			vector < vector <string> > dummy_vec;
			in_mapp_string.push_back(dummy_vec);
			in_mapp_string[tot_lig-1].push_back(in_record_tmp);
		}
	}

	//-> final collect
	int j,k;
	char command[5];
	int tot_num=0;
	for(i=0;i<tot_lig;i++)
	{
		for(j=0;j<in_mapp_int[i];j++)
		{
			//get name
			sprintf(command,"%-4d",tot_num);
			command[4]='\0';
			for(k=0;k<(int)strlen(command);k++)if(command[k]==' ')command[k]='_';
			string name_tmp=command;
			string name=in_mapp_nam[i]+"_"+name_tmp;
			//get type
			string wsss=in_mapp_nam[i];
			for(k=0;k<3;k++)if(wsss[k]=='_')wsss[k]=' ';
			int result=0;
			for(k=0;k<3;k++)result+=(Ligand_Trans(wsss[k]))*(int)pow(37.0,1.0*k);
			int retv=Ligand_Num_Code(result);
			const char *ligtype=Ligand_CAMEO(retv);
			//assign
			vector <XYZ> xyz_tmp;
			int moln=Ligand_String_To_XYZ(in_mapp_string[i][j],xyz_tmp);
			Ligand_Struc ligand_struc;
			ligand_struc.lig_name=name;
			ligand_struc.lig_type=ligtype;
			ligand_struc.lig_moln=moln;
			ligand_struc.lig_xyz=xyz_tmp;
			ligand_struc.lig_data=in_mapp_string[i][j];
			//incremental
			output.push_back(ligand_struc);
			tot_num++;
		}
	}
	//-> final return
	return tot_num;
}


//================== All_Chain_Process =============//__110530__//
//-> Compare ligand and chain 
int Compare_Ligand_and_Chain(vector <XYZ> &protein, vector <XYZ> &ligand,double r_cut)
{
	int i,j;
	int ll=(int)ligand.size();
	int pl=(int)protein.size();
	//distance check
	double dist2;
	double thres2=r_cut*r_cut;
	for(i=0;i<ll;i++)
	{
		for(j=0;j<pl;j++)
		{
			dist2=ligand[i].distance_square(protein[j]);
			if(dist2<thres2)return 1; //success !! (within radius)
		}
	}
	//final return
	return 0; //failed !! (not within radius)
}


//-> Compare ligand and chain complex
double Residue_Ligand_Distance(PDB_Residue & PDB, vector <XYZ> &ligand, double thres,int &atom_num)
{
	int i,k;
	int number;
	int totnum;
	XYZ xyz;
	double dist2;
	double minval=99999;
	totnum=(int)ligand.size();
	atom_num=0;
	double thres2=thres*thres;
	//backbone
	number=PDB.get_backbone_totnum();
	for(k=0;k<number;k++)
	{
		//get_backbone
		if(PDB.get_backbone_part_index(k)==0)continue;
		PDB.get_backbone_atom(k,xyz);
		//calculate distance
		for(i=0;i<totnum;i++)
		{
			dist2=ligand[i].distance_square(xyz);
			if(dist2<minval)minval=dist2;
			if(dist2<thres2)atom_num++;
		}
	}
	//sidechain
	number=PDB.get_sidechain_totnum();
	for(k=0;k<number;k++)
	{
		//get_sidechain
		if(PDB.get_sidechain_part_index(k)==0)continue;
		PDB.get_sidechain_atom(k,xyz);
		//calculate distance
		for(i=0;i<totnum;i++)
		{
			dist2=ligand[i].distance_square(xyz);
			if(dist2<minval)minval=dist2;
			if(dist2<thres2)atom_num++;
		}
	}
	//return
	return minval;
}

//--------- misc functions --------//start
void Kill_Space(string &in,string &out)
{
	int i;
	out="";
	for(i=0;i<(int)in.length();i++)
	{
		if(in[i]!=' ')out+=in[i];
	}
}
int Check_Ins(string &in)
{
	int i=(int)in.length()-1;
	if(in[i]>='0'&&in[i]<='9')return 0;
	else return 1;
}
int Check_Chain(string &in)
{
	int i=0;
	if(in[i]=='|')return 0;  //-> null chain
	else return 1;
}
//--------- misc functions --------//end

int Compare_Ligand_and_Chain_Complex(vector <PDB_Residue> &protein, vector <XYZ> &ligand,double r_cut, 
	vector <int> &pos_rec, vector <char> &cha_rec, vector <string> &ind_rec, vector <double> &min_rec, int &atom_num)
{
	int i;
	int pl=(int)protein.size();
	//distance check
	double dist2;
	double thres2=r_cut*r_cut;
	pos_rec.clear();
	cha_rec.clear();
	ind_rec.clear();
	min_rec.clear();
	int count=0;
	int cur_atom;
	atom_num=0;
	for(i=0;i<pl;i++)
	{
		dist2=Residue_Ligand_Distance(protein[i], ligand, r_cut, cur_atom);
		if(dist2<thres2)
		{
			pos_rec.push_back(i);
			char c=protein[i].get_AA();
			string ind_;
			protein[i].get_PDB_residue_number(ind_);
			string ind;
			ind=ind_.substr(1,5);
			Kill_Space(ind,ind_);
			cha_rec.push_back(c);
			ind_rec.push_back(ind_);
			min_rec.push_back(sqrt(dist2));
			count++;
			atom_num+=cur_atom;
		}
	}
	//return
	return count;
}


//================= for Point-Cloud ===================//__190602__//
//-> record point-cloud
void Residue_Ligand_Distance_PC(PDB_Residue &PDB, int posi_val, int pacc_val, vector <Ligand_Struc> &ligands, 
	double r_cut, vector <XYZ> &pc, vector <string> &atom, vector <int> &atom_lab,
	vector <string> &posi, vector <char> &resi, vector <char> &label, vector <int> &pacc, vector <double> &bfac)
{
	int i,j,k;
	int number;
	double dist2;
	double thres2=r_cut*r_cut;
	//b-factor
	XYZ xyz;
	int numb;
	double rfactor,temperature;
	//get amino
	char amino=PDB.get_AA();
	//get_backbone
	number=PDB.get_backbone_totnum();
	for(k=0;k<number;k++)
	{
		//get backbone atom
		if(PDB.get_backbone_part_index(k)==0)continue;
		XYZ xyz;
		PDB.get_backbone_atom(k,xyz);
		const char *atomname=backbone_atom_name_decode(k);
		string atom_name=atomname;
		string pdbind_;
		string posi_rel;
		PDB.get_PDB_residue_number(pdbind_);
		//-> reformat to A|123 style //-> start
		string chain="";
		chain.push_back(pdbind_[0]);
		string resi_prev=pdbind_.substr(1,pdbind_.length()-1);
		string posi_=chain+"|"+resi_prev;
		Kill_Space(posi_,posi_rel);
		if(Check_Ins(posi_rel)!=1)posi_rel.push_back(' ');
		if(Check_Chain(posi_rel)==0)posi_rel=" "+posi_rel;
		//-> reformat to A|123 style //-> end
		PDB.get_backbone_atom(k,xyz, numb, rfactor, temperature);
		//record as point-cloud
		pc.push_back(xyz);
		atom.push_back(atom_name);
		atom_lab.push_back(k);
		posi.push_back(posi_rel);
		resi.push_back(amino);
		pacc.push_back(pacc_val);
		bfac.push_back(temperature);
		//calculate distance
		int retv=0;
		char retl='X';
		for(i=0;i<(int)ligands.size();i++)
		{
			for(j=0;j<(int)ligands[i].lig_xyz.size();j++)
			{
				dist2=ligands[i].lig_xyz[j].distance_square(xyz);
				if(dist2<thres2)
				{
					retv=1;
					retl=ligands[i].lig_type[0];
					break;
				}
			}
		}
		if(retv==1)label.push_back(retl);
		else label.push_back('0');
	}
	//get_sidechain
	number=PDB.get_sidechain_totnum();
	for(k=0;k<number;k++)
	{
		//get backbone atom
		if(PDB.get_sidechain_part_index(k)==0)continue;
		XYZ xyz;
		PDB.get_sidechain_atom(k,xyz);
		const char *atomname=sidechain_atom_name_decode(k,amino);
		string atom_name=atomname;
		string pdbind_;
		string posi_rel;
		PDB.get_PDB_residue_number(pdbind_);
		//-> reformat to A|123 style //-> start
		string chain="";
		chain.push_back(pdbind_[0]);
		string resi_prev=pdbind_.substr(1,pdbind_.length()-1);
		string posi_=chain+"|"+resi_prev;
		Kill_Space(posi_,posi_rel);
		if(Check_Ins(posi_rel)!=1)posi_rel.push_back(' ');
		if(Check_Chain(posi_rel)==0)posi_rel=" "+posi_rel;
		//-> reformat to A|123 style //-> end
		PDB.get_sidechain_atom(k,xyz, numb, rfactor, temperature);
		//record as point-cloud
		pc.push_back(xyz);
		atom.push_back(atom_name);
		atom_lab.push_back(k+4);
		posi.push_back(posi_rel);
		resi.push_back(amino);
		pacc.push_back(pacc_val);
		bfac.push_back(temperature);
		//calculate distance
		int retv=0;
		char retl='X';
		for(i=0;i<(int)ligands.size();i++)
		{
			for(j=0;j<(int)ligands[i].lig_xyz.size();j++)
			{
				dist2=ligands[i].lig_xyz[j].distance_square(xyz);
				if(dist2<thres2)
				{
					retv=1;
					retl=ligands[i].lig_type[0];
					break;
				}
			}
		}
		if(retv==1)label.push_back(retl);
		else label.push_back('0');
	}
}

//-> extract ALL residues
void Residue_Ligand_PC(vector <PDB_Residue> &chain, vector <int> &pacc_val, vector <Ligand_Struc> &ligands, 
	double r_cut, vector <XYZ> &pc, vector <string> &atom, vector <int> &atom_lab,
	vector <string> &posi, vector <char> &resi, vector <char> &label, vector <int> &pacc, vector <double> &bfac)
{
	pc.clear();
	atom.clear();
	atom_lab.clear();
	posi.clear();
	resi.clear();
	label.clear();
	pacc.clear();
	bfac.clear();
	int number=(int)chain.size();
	for(int i=0;i<number;i++)
	{
		Residue_Ligand_Distance_PC(chain[i],i+1,pacc_val[i],ligands,r_cut,pc,atom,atom_lab,posi,resi,label,pacc,bfac);
	}
}

//-> calculate relative solvent accessibility
void Calculate_ACC(vector <PDB_Residue> &chain, vector <int> &pAcc)
{
	//get size
	int moln=(int)chain.size();
	pAcc.resize(moln);
	//create data_structure
	Acc_Surface *acc_surface=new Acc_Surface(moln);
	XYZ **mcc;
	NewArray2D(&mcc,moln,15);
	int *mcc_side=new int[moln];
	char *ami=new char[moln+1];
	//calculate required data_structure
	for(int i=0;i<moln;i++)
	{
		chain[i].get_XYZ_array(mcc[i],mcc_side[i]);
		char amino=chain[i].get_AA();
		ami[i]=amino;
	}
	ami[moln]='\0';
	//calculate acc and pAcc
	char *acc=new char[moln+1];
	acc_surface->AC_Calc_SolvAcc(mcc,ami,moln,acc,mcc_side);
	acc[moln]='\0';
	for(int i=0;i<moln;i++)pAcc[i]=acc_surface->AC_normal[i];
	//delete
	delete [] acc;
	delete [] ami;
	delete [] mcc_side;
	DeleteArray2D(&mcc,moln);
	delete acc_surface;
}


//-> filter point-cloud
void Filter_PointCloud(int type,
	vector <XYZ> &pc_in, vector <string> &atom_in, vector <int> &atomlab_in, vector <string> &posi_in, vector <char> &resi_in, vector <char> &label_in, vector <int> &pacc_in, vector <double> &bfac_in,
	vector <XYZ> &pc_out, vector <string> &atom_out, vector <int> &atomlab_out, vector <string> &posi_out, vector <char> &resi_out, vector <char> &label_out, vector <int> &pacc_out, vector <double> &bfac_out)
{
	pc_out.clear();
	atom_out.clear();
	atomlab_out.clear();
	posi_out.clear();
	resi_out.clear();
	label_out.clear();
	pacc_out.clear();
	bfac_out.clear();
	long number=(long)pc_in.size();
	for(long i=0;i<number;i++)
	{
		//get atom_name
		string atom_name=atom_in[i];
		//judge type
		if(type<0) //-> CA+CB
		{
			if(atom_name=="CA " || atom_name=="CB ")
			{
				pc_out.push_back(pc_in[i]);
				atom_out.push_back(atom_in[i]);
				atomlab_out.push_back(atomlab_in[i]);
				posi_out.push_back(posi_in[i]);
				resi_out.push_back(resi_in[i]);
				label_out.push_back(label_in[i]);
				pacc_out.push_back(pacc_in[i]);
				bfac_out.push_back(bfac_in[i]);
			}
		}
		else if(type==0) //-> backbone+CB
		{
			if(atom_name=="N  " || atom_name=="CA " || atom_name=="C  " || atom_name=="O  " || atom_name=="CB ")
			{
				pc_out.push_back(pc_in[i]);
				atom_out.push_back(atom_in[i]);
				atomlab_out.push_back(atomlab_in[i]);
				posi_out.push_back(posi_in[i]);
				resi_out.push_back(resi_in[i]);
				label_out.push_back(label_in[i]);
				pacc_out.push_back(pacc_in[i]);
				bfac_out.push_back(bfac_in[i]);
			}
		}
		else
		{
			//--- filter out CB in glycine ---//
			if(resi_in[i]=='G' && atom_name=="CB ")continue;
			pc_out.push_back(pc_in[i]);
			atom_out.push_back(atom_in[i]);
			atomlab_out.push_back(atomlab_in[i]);
			posi_out.push_back(posi_in[i]);
			resi_out.push_back(resi_in[i]);
			label_out.push_back(label_in[i]);
			pacc_out.push_back(pacc_in[i]);
			bfac_out.push_back(bfac_in[i]);
		}
	}
}

//-> output PC to file
void Output_File_PC(FILE *fp,
	vector <XYZ> &pc, vector <string> &atom, vector <int> &atom_lab,
	vector <string> &posi, vector <char> &resi, vector <char> &label, vector <int> &pacc, vector <double> &bfac)
{
	for(long i=0;i<(long)pc.size();i++)
	{
		fprintf(fp,"%7s %c%c%c %8.3f %8.3f %8.3f %c %5.2f %3d\n",
			posi[i].c_str(),resi[i],atom[i][0],atom_lab[i]+'a',pc[i].X,pc[i].Y,pc[i].Z,label[i],bfac[i],pacc[i]);
	}
}


//---------------------- reconstruct missing heavy-atom --------------------------//
int Reconstruct_Missing(vector <PDB_Residue> &pdb)
{
	//init check
	int i,j;
	int moln=(int)pdb.size();
	for(i=0;i<moln;i++)if(pdb[i].get_backbone_part_index(1)==0)return -1; //no CA !!
	//data
	Confo_Beta confo_beta(moln);
	Confo_Back confo_back(moln);
	Confo_Lett confo_lett;
	XYZ *mol=new XYZ[moln];      //CA
	char *ami=new char[moln+1];  //ami
	//assign
	char amino;
	for(i=0;i<moln;i++)
	{
		pdb[i].get_backbone_atom(1,mol[i]);
		amino=pdb[i].get_AA();
		ami[i]=amino;
	}
	ami[moln]='\0';
	//check
	int correct;
	int iret;
	correct=1; //default:OK
	for(i=0;i<moln;i++)
	{
		iret=pdb[i].PDB_residue_backbone_check(4);
		if(iret!=1)
		{
			correct=0;
			break;
		}
		iret=pdb[i].PDB_residue_CB_check();
		if(iret!=1)
		{
			correct=0;
			break;
		}
	}
	if(correct==1)
	{
		delete [] mol;
		delete [] ami;
		return 0; //no modification
	}
	//reconstruct
	//[0]data
	char *cle=new char[moln+1];
	XYZ *mcb=new XYZ[moln];     //CB
	XYZ **mbb;                  //BackBone (N,CA,C,O,CB)
	NewArray2D(&mbb,moln,5);
	//[1]recon
	confo_lett.btb_ori(0,0,0,moln,mol,cle);
	cle[moln]='\0';
	confo_back.Recon_Back_Main(mol,cle,moln,mbb);      //given CA, recon BackBone (N,CA,C,O,CB)
	confo_beta.Recon_Beta_21(mol,mcb,moln,ami,cle);    //given CA, recon CB
	//[2]assign
	for(i=0;i<moln;i++)
	{
		//backbone (N,CA,C,O)
		for(j=0;j<4;j++)
		{
			if(pdb[i].get_backbone_part_index(j)==0)
			{
				pdb[i].set_backbone_atom(j,mbb[i][j]);
			}
		}
		//sidechain (CB)
		if(pdb[i].get_sidechain_part_index(0)==0)
		{
			if(i==0||i==moln-1)pdb[i].set_sidechain_atom(0,mbb[i][4]);
			else pdb[i].set_sidechain_atom(0,mcb[i]);
		}
	}
	//final return
	delete [] mol;
	delete [] ami;
	delete [] mcb;
	delete [] cle;
	DeleteArray2D(&mbb,moln);
	return 1; //reconstruct
}

//-> PDB_Residue_To_XYZ
int PDB_Residue_To_XYZ(PDB_Residue &PDB, vector <XYZ> &output)
{
	int k;
	int number;
	XYZ xyz;
	int count=0;
	//backbone_out
	number=PDB.get_backbone_totnum();  //this value should be 4
	for(k=0;k<number;k++)
	{
		//get_backbone
		if(PDB.get_backbone_part_index(k)==0)continue;
		PDB.get_backbone_atom(k,xyz);
		output.push_back(xyz);
		count++;
	}
	char amino=PDB.get_AA();
	if(amino=='G')return count;
	//sidechain_out
	number=PDB.get_sidechain_totnum();
	for(k=0;k<number;k++)
	{
		//get_sidechain
		if(PDB.get_sidechain_part_index(k)==0)continue;
		PDB.get_sidechain_atom(k,xyz);
		output.push_back(xyz);
		count++;
	}
	//return
	return count;
}
int PDB_Residue_To_XYZ_Total(vector <PDB_Residue> &pdb, vector <XYZ> &output)
{
	int i;
	int count=0;
	int retv=0;
	int moln=(int)pdb.size();
	output.clear();
	for(i=0;i<moln;i++)
	{
		retv=PDB_Residue_To_XYZ(pdb[i],output);
		count+=retv;
	}
	return count;
}

//--------- filter ligands -----------//
//-> example
/*
HOH
SO4
...
*/
void Filter_Ligands(string &infile, 
	vector <Ligand_Struc> &ligands_in,
	vector <Ligand_Struc> &ligands_out)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(infile.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",infile.c_str());
		exit(-1);
	}
	//proc
	vector <string> filter_list;
	map <string, int > name_mapping;
	map<string, int >::iterator iter;
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		if(buf=="")continue;
		if(buf.length()<3)continue;
		temp=buf.substr(0,3);
		iter = name_mapping.find(temp);
		if(iter == name_mapping.end())
		{
			count++;
			name_mapping.insert(map < string, int >::value_type(temp, count));
			filter_list.push_back(temp);
		}
	}
	//calc
	ligands_out.clear();
	for(int i=0;i<(int)ligands_in.size();i++)
	{
		string lig_name=ligands_in[i].lig_name.substr(0,3);
		iter = name_mapping.find(lig_name);
		if(iter == name_mapping.end())
		{
			ligands_out.push_back(ligands_in[i]);
		}
	}
}


//-> Extract all chains for ligands
//[note]: len_cut -> for pdb_chain
//        r_cut   -> for ligand
int PDB_Ligand_All_Process(string &file,string &out_name,
	string &pdb_out_dir,string &ligand_out_dir,string &log_out_dir,string &pc_out_dir,
	int len_cut,double r_cut)
{
	//class
	Mol_File mol_input;
	Mol_Out mol_output;
	mol_input.MODRES=1;
	//data
	vector <PDB_Residue> pdb;
	vector <XYZ> mol;
	vector <XYZ> total_xyz;
	string ami;
	vector <string> ind;
	PDB_Chain pdb_chain;
	PDB_Residue residue;
	int retv;
	//chain load
	vector <PDB_Chain> chains;
	retv=mol_input.PDB_read_pdb_file(file,chains);
	if(retv<0)return 0; //failed
	//ligand load	
	vector <Ligand_Struc> ligands;
	retv=PDB_Extract_Ligand(file,ligands);
	if(retv<0)return 0; //failed

	//======================== filter ligands according to instructions ======================//start
	//--- filter-out ligands ----//type
	if(LG_TYPE!="")
	{
		vector <Ligand_Struc> ligands_filter;
		for(int i=0;i<(int)ligands.size();i++)
		{
			char lig_type=ligands[i].lig_type[0];
			for(int k=0;k<(int)LG_TYPE.length();k++)
			{
				if(LG_TYPE[k]==lig_type)
				{
					ligands_filter.push_back(ligands[i]);
					break;
				}
			}
		}
		ligands=ligands_filter;
	}
	//--- filter-out ligands ----//name
	if(LG_FILTER!="")
	{
		vector <Ligand_Struc> ligands_filter;
		Filter_Ligands(LG_FILTER,ligands,ligands_filter);
		ligands=ligands_filter;
	}
	//--- filter-out ligands ----//number
	if(LG_MINNUM>1)
	{
		vector <Ligand_Struc> ligands_filter;
		for(int i=0;i<(int)ligands.size();i++)
		{
			if(ligands[i].lig_moln < LG_MINNUM)continue;
			ligands_filter.push_back(ligands[i]);
		}
		ligands=ligands_filter;
	}
	//======================== filter ligands according to instructions ======================//end

	//--- output ligand ----//
	{
		for(int i=0;i<(int)ligands.size();i++)
		{
			string cur_nam=out_name+"_"+ligands[i].lig_name;
			string outname=ligand_out_dir+"/"+cur_nam+".pdb";
			FILE *fp=fopen(outname.c_str(),"wb");
			for(int k=0;k<(int)ligands[i].lig_data.size();k++)
			{
				fprintf(fp,"%s\n",ligands[i].lig_data[k].c_str());
			}
			fclose(fp);
		}
		if(LOG_OR_NOT==1 && (int)ligands.size()>0)
		{
			string log_file=log_out_dir+"/"+out_name+".ligand_size";
			FILE *f2=fopen(log_file.c_str(),"wb");
			for(int i=0;i<(int)ligands.size();i++)
			{
				string cur_nam=out_name+"_"+ligands[i].lig_name;
				fprintf(f2,"%s %d\n",cur_nam.c_str(),ligands[i].lig_moln);
			}
			fclose(f2);
		}
	}

	//chain-ligand process
	int i,j,k;
	int moln;
	string chain;
	string cur_nam;
	string resnum;
	int chain_size;
	chain_size=(int)chains.size();
	int lig_size;
	lig_size=(int)ligands.size();
	//fprintf data
	FILE *fp;
	int has_ligand;
	string outname;
	//log
	FILE *f1,*f2,*f3,*f4,*f5;
	string log_file;
	f1=0;
	f2=0;
	f3=0;
	f4=0;
	f5=0;
	//real process
	int ligand_log_first=1;
	int chain_log_first=1;
	for(i=0;i<chain_size;i++)
	{
		//read
		pdb_chain=chains[i];
		moln=pdb_chain.get_length();
		chain=pdb_chain.get_chain_id();
		if(chain==" ")chain="_";
		//length check
		if(moln<len_cut)continue;
		//assign
		pdb.resize(moln);
		mol.resize(moln);
		ami.resize(moln);
		ind.resize(moln);
		for(k=0;k<moln;k++)
		{
			pdb_chain.get_residue(k,residue);
			pdb[k]=residue;
			residue.get_backbone_atom(1,mol[k]);
			ami[k]=residue.get_AA();
			residue.get_PDB_residue_number(resnum);
			resnum=resnum.substr(1,5);
			ind[k]=resnum;
		}
		Reconstruct_Missing(pdb);
		PDB_Residue_To_XYZ_Total(pdb,total_xyz); 
		//check length and ligand
		has_ligand=0;
		for(j=0;j<lig_size;j++)
		{
			retv=Compare_Ligand_and_Chain(total_xyz,ligands[j].lig_xyz,r_cut);
			if(retv==0)continue;
			//get detailed binding site
			vector <int> pos_rec;
			vector <char> cha_rec;
			vector <string> ind_rec;
			vector <double> min_rec;
			int atom_count=0;
			int count=Compare_Ligand_and_Chain_Complex(pdb,ligands[j].lig_xyz,r_cut,pos_rec,cha_rec,ind_rec,min_rec,atom_count);
			//output ligand_log
			if(LOG_OR_NOT==1 && count>=LG_MINCUT && atom_count>=LG_ATOMCUT)
			{
				if(ligand_log_first==1)
				{
					ligand_log_first=0;
					log_file=log_out_dir+"/"+out_name+".ligand_log";
					f2=fopen(log_file.c_str(),"wb");
					fprintf(f2,">%s\n",out_name.c_str());
					log_file=log_out_dir+"/"+out_name+".ligand_chain";
					f3=fopen(log_file.c_str(),"wb");
				}
				cur_nam=out_name+chain+"_"+ligands[j].lig_name;
				fprintf(f2,"%s -> ",cur_nam.c_str());
				for(k=0;k<count;k++)fprintf(f2,"%d|%s|%c|%3.1f ",pos_rec[k]+1,ind_rec[k].c_str(),cha_rec[k],min_rec[k]);
				fprintf(f2,"\n");
				//has ligand
				has_ligand++;
			}
		}
		if(LOG_OR_NOT==1)
		{
			if(has_ligand>0)fprintf(f3,"%s%s\n",out_name.c_str(),chain.c_str());
		}
		//output pdb
		cur_nam=out_name+chain;
		outname=pdb_out_dir+"/"+cur_nam+".pdb";
		fp=fopen(outname.c_str(),"wb");
		mol_output.Output_PDB_III(fp,moln,pdb,chain[0]);
		fclose(fp);
		//output pdb_log
		if(LOG_OR_NOT==1)
		{
			if(chain_log_first==1)
			{
				chain_log_first=0;
				log_file=log_out_dir+"/"+out_name+".chain_log";
				f1=fopen(log_file.c_str(),"wb");
				log_file=log_out_dir+"/"+out_name+".atom_seq";
				f4=fopen(log_file.c_str(),"wb");
				log_file=log_out_dir+"/"+out_name+".atom_ind";
				f5=fopen(log_file.c_str(),"wb");
			}
			fprintf(f1,"%s %4d\n",cur_nam.c_str(),moln);
			fprintf(f4,">%s\n",cur_nam.c_str());
			fprintf(f4,"%s\n",ami.c_str());
			fprintf(f5,">%s\n",cur_nam.c_str());
			for(k=0;k<moln;k++)fprintf(f5,"%s",ind[k].c_str());
			fprintf(f5,"\n");
		}
	}
	//final
	if(LOG_OR_NOT==1)
	{
		if(ligand_log_first==0)
		{
			fclose(f2);
			fclose(f3);
		}
		if(chain_log_first==0)
		{
			fclose(f1);
			fclose(f4);
			fclose(f5);
		}
	}

	//================= record binding residues in Point-Cloud ==================//__190602__//
	if(pc_out_dir!="")
	{
		//--- main process ---//
		for(i=0;i<chain_size;i++)
		{
			//read
			pdb_chain=chains[i];
			moln=pdb_chain.get_length();
			chain=pdb_chain.get_chain_id();
			if(chain==" ")chain="_";
			//length check
			if(moln<len_cut)continue;
			//assign
			pdb.resize(moln);
			for(k=0;k<moln;k++)
			{
				pdb_chain.get_residue(k,residue);
				pdb[k]=residue;
			}
			Reconstruct_Missing(pdb);
			//calculate pAcc
			vector <int> pacc;
			Calculate_ACC(pdb,pacc);
			//calculate point-cloud
			vector <XYZ> pc_in;
			vector <string> atom_in;
			vector <int> atomlab_in;
			vector <string> posi_in;
			vector <char> resi_in;
			vector <char> label_in;
			vector <int> pacc_in;
			vector <double> bfac_in;
			Residue_Ligand_PC(pdb,pacc,ligands,r_cut,pc_in,atom_in,atomlab_in,
				posi_in,resi_in,label_in,pacc_in,bfac_in);
			//filter point-cloud according to type
			vector <XYZ> pc_out;
			vector <string> atom_out;
			vector <int> atomlab_out;
			vector <string> posi_out;
			vector <char> resi_out;
			vector <char> label_out;
			vector <int> pacc_out;
			vector <double> bfac_out;
			Filter_PointCloud(PC_TYPE,
				pc_in,atom_in,atomlab_in,posi_in,resi_in,label_in,pacc_in,bfac_in,
				pc_out,atom_out,atomlab_out,posi_out,resi_out,label_out,pacc_out,bfac_out);
			//output point-cloud to file
			cur_nam=out_name+chain; 
			outname=pc_out_dir+"/"+cur_nam+".pc_xyz";
			fp=fopen(outname.c_str(),"wb");
			Output_File_PC(fp,pc_out,atom_out,atomlab_out,posi_out,resi_out,label_out,pacc_out,bfac_out);
			fclose(fp);
		}
	}

	return 1; //success
}

//----- check ligand type string ----//
int Check_LGtype(string &in)
{
	for(int i=0;i<(int)in.length();i++)
	{
		if(Ligand_Type(in[i])<0)
		{
			fprintf(stderr,"-t LG_TYPE %s should contain 'I', 'O', 'P', 'N', or 'X' \n",in.c_str());
			exit(-1);
		}
	}
}


//============== main ===============//
int main(int argc, char** argv)
{
	//---- PDB_Ligand_Process ---//process all chains and ligands
	{
		process_args(argc,argv);
		//judge
		if(INPUT_FILE=="")
		{
			fprintf(stderr,"-i INPUT_FILE should not be blank \n");
			exit(-1);
		}
		if(OUTPUT_NAME=="")
		{
			getBaseName(INPUT_FILE,OUTPUT_NAME,'/','.');
		}
		if(LENGTH_CUT<0)
		{
			fprintf(stderr,"-n LENGTH_CUT %d should be positive \n",LENGTH_CUT);
			exit(-1);
		}
		if(DISTANCE_CUT<0)
		{
			fprintf(stderr,"-d DISTANCE_CUT %lf should be positive \n",DISTANCE_CUT);
			exit(-1);
		}
		if(LOG_OR_NOT<0 || LOG_OR_NOT>1)
		{
			fprintf(stderr,"-l LOG_OR_NOT %d should be 0 or 1 \n",LOG_OR_NOT);
			exit(-1);
		}
		if(PC_TYPE<-1 || PC_TYPE>1)
		{
			fprintf(stderr,"-T PC_TYPE %d should be -1,0,1 \n",PC_TYPE);
			exit(-1);
		}
		Check_LGtype(LG_TYPE);
		//assign
		string input_pdb=INPUT_FILE;
		string out_name=OUTPUT_NAME;
		string pdb_out_root=CHAIN_OUTROOT;
		string ligand_out_root=LIGAND_OUTROOT;
		string log_out_root=LOG_ROOT;
		string point_cloud_root=PC_ROOT;
		//debug
		if(LOG_OR_NOT==1)
		{
			if(log_out_root=="")log_out_root="./";
		}
		//create output
		char command[30000];
		sprintf(command,"mkdir -p %s",pdb_out_root.c_str());
		system(command);
		sprintf(command,"mkdir -p %s",ligand_out_root.c_str());
		system(command);
		if(log_out_root!="")
		{
			sprintf(command,"mkdir -p %s",log_out_root.c_str());
			system(command);
		}
		if(point_cloud_root!="")
		{
			sprintf(command,"mkdir -p %s",point_cloud_root.c_str());
			system(command);
		}
		//process
		int len_cut=LENGTH_CUT;
		double distance_cut=DISTANCE_CUT;
		PDB_Ligand_All_Process(input_pdb,out_name,
			pdb_out_root,ligand_out_root,log_out_root,point_cloud_root,
			len_cut,distance_cut);
		exit(0);
	}
}




