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

//--------- Parse_Str_Str --------//
int Parse_Str_Str(string &in,vector <string> &out, char separator)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(!getline(www,buf,separator))break;
		out.push_back(buf);
		count++;
	}
	return count;
}
//--------- Parse_Str_Str (automatic) --------//
int Parse_Str_Str(string &in,vector <string> &out)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(! (www>>buf) )break;
		out.push_back(buf);
		count++;
	}
	return count;
}

//------- int to string --------//
template <typename T>
string NumberToString( T Number )
{
	ostringstream ss;
	ss << Number;
	return ss.str();
}



//================================== FEATURE GENERATION FOR XYZ FILE =============================//
//-> We use the atom-level and residue-level features according to the following paper:
// P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure

//=============== atom features =================//
const int Resi_Number=20;
const int Atom_Number=167;
double Atom_Features[Atom_Number][13]={
//selection:    1             1            0              0                 1            0   0        1        1        1          1           1         1
//atomName,apRawValids,apRawInvalids,ap5sasaValids,ap5sasaInvalids,atomicHydrophobicity,aph,anh, vsAromatic,vsCation,vsAnion,vsHydrophobic,vsAcceptor,vsDonor
{/*ALA.N,   */ 0.772020725,0.617283951,0.851762434,0.513580247,1,1,0                          , 0,0,0,0,0,1 },
{/*ALA.CA,  */ 0.627277286,0.514029181,0.744567842,0.555555556,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ALA.C,   */ 0.414674912,0.351290685,0.442781265,0.372839506,1,1,0                          , 0,0,0,1,0,0 },
{/*ALA.O,   */ 0.7100117,0.554433221,0.785127958,0.567901235,1,1,0                            , 0,0,0,0,1,0 },
{/*ALA.CB,  */ 1.551729901,1.052749719,2.233703525,1.241975309,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ARG.N,   */ 0.320580674,0.210679612,0.376160991,0.219217491,1,1,0                          , 0,0,0,0,0,1 },
{/*ARG.CA,  */ 0.234840466,0.166990291,0.269543344,0.181818182,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ARG.C,   */ 0.157568426,0.115533981,0.175696594,0.116225547,1,1,0                          , 0,0,0,1,0,0 },
{/*ARG.O,   */ 0.202631181,0.177184466,0.229295666,0.186996548,1,1,0                          , 0,0,0,0,1,0 },
{/*ARG.CB,  */ 0.380311508,0.263592233,0.449690402,0.285960875,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ARG.CG,  */ 0.425979132,0.276699029,0.507352941,0.306098964,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ARG.CD,  */ 0.708150612,0.430097087,0.816756966,0.466052934,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ARG.NE,  */ 0.988658703,0.697572816,1.099845201,0.722669735,1,1,0                          , 0,1,0,0,0,0 },
{/*ARG.CZ,  */ 1.332678058,0.988349515,1.353134675,0.98273878,1,1,0                           , 0,0,0,1,0,0 },
{/*ARG.NH1, */ 1.6215031,1.188349515,1.621904025,1.180667434,1,1,0                            , 0,1,0,0,0,0 },
{/*ARG.NH2, */ 2.04929684,1.574757282,2.061532508,1.54027618,1,1,0                            , 0,1,0,0,0,0 },
{/*ASN.N,   */ 0.422993492,0.323799796,0.589206066,0.37,1,1,0                                 , 0,0,0,0,0,1 },
{/*ASN.CA,  */ 0.327765727,0.31052094,0.411685995,0.35,-1,0,-1                                , 0,0,0,1,0,0 },
{/*ASN.C,   */ 0.278524946,0.251276813,0.328724353,0.252857143,1,1,0                          , 0,0,0,1,0,0 },
{/*ASN.O,   */ 0.45032538,0.335035751,0.433987511,0.372857143,1,1,0                           , 0,0,0,0,1,0 },
{/*ASN.CB,  */ 0.536442516,0.592441267,0.712756467,0.664285714,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ASN.CG,  */ 0.934707158,0.93258427,1.047725245,1.092857143,1,1,0                           , 0,0,0,1,0,0 },
{/*ASN.OD1, */ 1.353145336,0.916241062,1.475022302,1.044285714,1,1,0                          , 0,0,0,0,1,0 },
{/*ASN.ND2, */ 1.890672451,1.816138917,2.076271186,2.018571429,1,1,0                          , 0,0,0,0,0,1 },
{/*ASP.N,   */ 0.334101382,0.358892439,0.478428622,0.499245852,1,1,0                          , 0,0,0,0,0,1 },
{/*ASP.CA,  */ 0.392059553,0.349307774,0.524026657,0.401206637,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ASP.C,   */ 0.270294222,0.261980831,0.37109786,0.286576169,1,1,0                           , 0,0,0,1,0,0 },
{/*ASP.O,   */ 0.290677065,0.38658147,0.407927043,0.408748115,1,1,0                           , 0,0,0,0,1,0 },
{/*ASP.CB,  */ 0.383551932,0.469648562,0.535250789,0.58974359,-1,0,-1                         , 0,0,0,1,0,0 },
{/*ASP.CG,  */ 1.321517192,0.774227902,1.27148369,0.690799397,1,1,0                           , 0,0,0,1,0,0 },
{/*ASP.OD1, */ 1.92786246,1.030883919,1.905296387,0.915535445,1,1,0                           , 0,0,1,0,0,0 },
{/*ASP.OD2, */ 2.158808933,1.020234292,1.923886356,0.859728507,1,1,0                          , 0,0,1,0,0,0 },
{/*CYS.N,   */ 0.492625369,0.39084507,0.688461538,0.5,1,1,0                                   , 0,0,0,0,0,1 },
{/*CYS.CA,  */ 0.597640118,0.788732394,0.730769231,1.082278481,-1,0,-1                        , 0,0,0,1,0,0 },
{/*CYS.C,   */ 0.464306785,0.221830986,0.392307692,0.151898734,1,1,0                          , 0,0,0,1,0,0 },
{/*CYS.O,   */ 0.791150442,0.464788732,0.507692308,0.17721519,1,1,0                           , 0,0,0,0,1,0 },
{/*CYS.CB,  */ 0.895575221,1.559859155,1.682692308,2.272151899,-1,0,-1                        , 0,0,0,1,0,0 },
{/*CYS.SG,  */ 1.424778761,2.637323944,2.478846154,3.759493671,1,1,0                          , 0,0,0,0,1,0 },
{/*GLN.N,   */ 0.314408911,0.355987055,0.405134627,0.413865546,1,1,0                          , 0,0,0,0,0,1 },
{/*GLN.CA,  */ 0.312971613,0.292880259,0.420162805,0.306722689,-1,0,-1                        , 0,0,0,1,0,0 },
{/*GLN.C,   */ 0.233201581,0.182847896,0.281778334,0.180672269,1,1,0                          , 0,0,0,1,0,0 },
{/*GLN.O,   */ 0.3729788,0.323624595,0.445835942,0.334033613,1,1,0                            , 0,0,0,0,1,0 },
{/*GLN.CB,  */ 0.438016529,0.454692557,0.60989355,0.474789916,-1,0,-1                         , 0,0,0,1,0,0 },
{/*GLN.CG,  */ 0.561264822,0.498381877,0.775829681,0.527310924,-1,0,-1                        , 0,0,0,1,0,0 },
{/*GLN.CD,  */ 0.982752425,0.558252427,1.065122104,0.56092437,1,1,0                           , 0,0,0,1,0,0 },
{/*GLN.OE1, */ 1.518864535,0.665048544,1.476518472,0.682773109,1,1,0                          , 0,0,0,0,1,0 },
{/*GLN.NE2, */ 1.584980237,1.116504854,1.57294928,1.031512605,1,1,0                           , 0,0,0,0,0,1 },
{/*GLU.N,   */ 0.208234737,0.288844622,0.306601467,0.405953992,1,1,0                          , 0,0,0,0,0,1 },
{/*GLU.CA,  */ 0.244202556,0.300796813,0.297310513,0.370771313,-1,0,-1                        , 0,0,0,1,0,0 },
{/*GLU.C,   */ 0.238050166,0.192231076,0.297310513,0.219215156,1,1,0                          , 0,0,0,1,0,0 },
{/*GLU.O,   */ 0.379081874,0.309760956,0.542787286,0.296346414,1,1,0                          , 0,0,0,0,1,0 },
{/*GLU.CB,  */ 0.214150497,0.330677291,0.320293399,0.414073072,-1,0,-1                        , 0,0,0,1,0,0 },
{/*GLU.CG,  */ 0.337908187,0.342629482,0.503178484,0.397834912,-1,0,-1                        , 0,0,0,1,0,0 },
{/*GLU.CD,  */ 1.122574539,0.707171315,1.012224939,0.638700947,1,1,0                          , 0,0,0,1,0,0 },
{/*GLU.OE1, */ 1.715097018,0.958167331,1.435696822,0.810554804,1,1,0                          , 0,0,1,0,0,0 },
{/*GLU.OE2, */ 1.835068623,1.070717131,1.561369193,0.964817321,1,1,0                          , 0,0,1,0,0,0 },
{/*GLY.N,   */ 1.417251755,0.84040404,1.717056475,1.208474576,1,1,0                           , 0,0,0,0,0,1 },
{/*GLY.CA,  */ 1.974833592,1.462626263,2.376497433,1.437288136,-1,0,-1                        , 0,0,0,1,0,0 },
{/*GLY.C,   */ 0.967630163,0.781144781,1.123217342,0.694915254,1,1,0                          , 0,0,0,1,0,0 },
{/*GLY.O,   */ 0.928877542,0.916498316,1.26982316,0.738983051,1,1,0                           , 0,0,0,0,1,0 },
{/*HIS.N,   */ 0.200358986,0.182815356,0.275819672,0.243172952,1,1,0                          , 0,0,0,0,0,1 },
{/*HIS.CA,  */ 0.229302221,0.31535649,0.318442623,0.358907672,-1,0,-1                         , 0,0,0,1,0,0 },
{/*HIS.C,   */ 0.184428988,0.18738574,0.246311475,0.204161248,1,1,0                           , 0,0,0,1,0,0 },
{/*HIS.O,   */ 0.309625309,0.219378428,0.377459016,0.183355007,1,1,0                          , 0,0,0,0,1,0 },
{/*HIS.CB,  */ 0.434821629,0.438756856,0.647131148,0.485045514,-1,0,-1                        , 0,0,0,1,0,0 },
{/*HIS.CG,  */ 0.442001346,0.363802559,0.639344262,0.425227568,-1,0,-1                        , 1,0,0,0,0,0 },
{/*HIS.ND1, */ 0.855732556,0.805301645,1.170901639,0.945383615,1,1,0                          , 0,0,0,0,1,1 },
{/*HIS.CD2, */ 0.981826341,1.103290676,1.158606557,1.301690507,-1,0,-1                        , 1,0,0,0,0,0 },
{/*HIS.CE1, */ 1.44940543,1.510968921,1.722131148,1.77893368,1,1,0                            , 1,0,0,0,0,0 },
{/*HIS.NE2, */ 1.792910029,2.045703839,1.961885246,2.356306892,1,1,0                          , 0,0,0,0,1,1 },
{/*ILE.N,   */ 0.450892857,0.393113343,0.47948244,0.368983957,1,1,0                           , 0,0,0,0,0,1 },
{/*ILE.CA,  */ 0.313108766,0.24964132,0.325693161,0.229946524,-1,0,-1                         , 0,0,0,1,0,0 },
{/*ILE.C,   */ 0.243912338,0.232424677,0.182994455,0.187165775,1,1,0                          , 0,0,0,1,0,0 },
{/*ILE.O,   */ 0.553977273,0.424677188,0.384473198,0.363636364,1,1,0                          , 0,0,0,0,1,0 },
{/*ILE.CB,  */ 0.343141234,0.241032999,0.436229205,0.278074866,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ILE.CG1, */ 0.669642857,0.304160689,0.931608133,0.409090909,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ILE.CG2, */ 0.974837662,0.522238164,1.241774492,0.582887701,-1,0,-1                        , 0,0,0,1,0,0 },
{/*ILE.CD1, */ 1.122767857,0.680057389,1.594085028,0.818181818,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LEU.N,   */ 0.476220672,0.470860315,0.489297475,0.274143302,1,1,0                          , 0,0,0,0,0,1 },
{/*LEU.CA,  */ 0.289949271,0.212765957,0.296103183,0.205607477,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LEU.C,   */ 0.305960685,0.222016651,0.273874863,0.177570093,1,1,0                          , 0,0,0,1,0,0 },
{/*LEU.O,   */ 0.593532023,0.37280296,0.459659715,0.387850467,1,1,0                           , 0,0,0,0,1,0 },
{/*LEU.CB,  */ 0.417564997,0.314523589,0.472557629,0.306853583,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LEU.CG,  */ 0.309923906,0.130434783,0.394072448,0.152647975,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LEU.CD1, */ 0.99445149,0.548566142,1.326289791,0.691588785,-1,0,-1                         , 0,0,0,1,0,0 },
{/*LEU.CD2, */ 0.887761573,0.624421832,1.107299671,0.721183801,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LYS.N,   */ 0.639755767,0.285328533,0.536626344,0.323809524,1,1,0                          , 0,0,0,0,0,1 },
{/*LYS.CA,  */ 0.397331524,0.246624662,0.442540323,0.247619048,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LYS.C,   */ 0.206919946,0.177317732,0.215053763,0.175661376,1,1,0                          , 0,0,0,1,0,0 },
{/*LYS.O,   */ 0.217548621,0.272727273,0.290994624,0.232804233,1,1,0                          , 0,0,0,0,1,0 },
{/*LYS.CB,  */ 0.44843962,0.334833483,0.487567204,0.364021164,-1,0,-1                         , 0,0,0,1,0,0 },
{/*LYS.CG,  */ 0.490275893,0.246624662,0.622647849,0.24973545,-1,0,-1                         , 0,0,0,1,0,0 },
{/*LYS.CD,  */ 0.557892356,0.508550855,0.704637097,0.532275132,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LYS.CE,  */ 1.342152872,0.884788479,1.356854839,0.814814815,-1,0,-1                        , 0,0,0,1,0,0 },
{/*LYS.NZ,  */ 2.42152872,1.596759676,2.103494624,1.453968254,1,1,0                           , 0,1,0,0,0,0 },
{/*MET.N,   */ 0.503381234,0.386075949,0.556343019,0.491071429,1,1,0                          , 0,0,0,0,0,1 },
{/*MET.CA,  */ 0.373203719,0.341772152,0.397590361,0.401785714,-1,0,-1                        , 0,0,0,1,0,0 },
{/*MET.C,   */ 0.256551141,0.246835443,0.265768958,0.272321429,1,1,0                          , 0,0,0,1,0,0 },
{/*MET.O,   */ 0.469146238,0.363924051,0.439404678,0.258928571,1,1,0                          , 0,0,0,0,1,0 },
{/*MET.CB,  */ 0.622992392,0.525316456,0.747696669,0.598214286,-1,0,-1                        , 0,0,0,1,0,0 },
{/*MET.CG,  */ 0.750211327,0.471518987,0.931963147,0.491071429,-1,0,-1                        , 0,0,0,1,0,0 },
{/*MET.SD,  */ 0.945477599,0.575949367,1.151665485,0.696428571,1,1,0                          , 0,0,0,1,0,0 },
{/*MET.CE,  */ 1.154691462,0.724683544,1.442948264,0.785714286,-1,0,-1                        , 0,0,0,1,0,0 },
{/*PHE.N,   */ 0.28952689,0.407496977,0.292449923,0.314102564,1,1,0                           , 0,0,0,0,0,1 },
{/*PHE.CA,  */ 0.223210675,0.296251511,0.245608629,0.271367521,-1,0,-1                        , 0,0,0,1,0,0 },
{/*PHE.C,   */ 0.215932066,0.23458283,0.212326656,0.226495726,1,1,0                           , 0,0,0,1,0,0 },
{/*PHE.O,   */ 0.432066316,0.336154776,0.388289676,0.39957265,1,1,0                           , 0,0,0,0,1,0 },
{/*PHE.CB,  */ 0.497573797,0.426844015,0.618181818,0.474358974,-1,0,-1                        , 0,0,0,1,0,0 },
{/*PHE.CG,  */ 0.45430651,0.209189843,0.62742681,0.277777778,-1,0,-1                          , 1,0,0,0,0,0 },
{/*PHE.CD1, */ 0.718358269,0.372430472,0.955932203,0.542735043,-1,0,-1                        , 1,0,0,0,0,0 },
{/*PHE.CD2, */ 0.709664375,0.340991536,0.930354391,0.405982906,-1,0,-1                        , 1,0,0,0,0,0 },
{/*PHE.CE1, */ 1.012939749,0.53808948,1.30385208,0.784188034,-1,0,-1                          , 1,0,0,0,0,0 },
{/*PHE.CE2, */ 1.000202184,0.498186215,1.2394453,0.566239316,-1,0,-1                          , 1,0,0,0,0,0 },
{/*PHE.CZ,  */ 1.160129397,0.529625151,1.422804314,0.698717949,-1,0,-1                        , 1,0,0,0,0,0 },
{/*PRO.N,   */ 0.235654998,0.413447783,0.277676951,0.220675944,-1,0,-1                        , 0,0,0,0,0,1 },
{/*PRO.CA,  */ 0.359076146,0.427753934,0.415003025,0.445328032,-1,0,-1                        , 0,0,0,1,0,0 },
{/*PRO.C,   */ 0.37531577,0.360515021,0.325468845,0.349900596,1,1,0                           , 0,0,0,1,0,0 },
{/*PRO.O,   */ 0.76290148,0.230329041,0.53599516,0.423459245,1,1,0                            , 0,0,0,0,1,0 },
{/*PRO.CB,  */ 0.628293035,0.519313305,0.698124622,0.562624254,-1,0,-1                        , 0,0,0,1,0,0 },
{/*PRO.CG,  */ 1.007578492,0.56223176,1.219600726,0.566600398,-1,0,-1                         , 0,0,0,1,0,0 },
{/*PRO.CD,  */ 0.985564778,0.555078684,1.199637024,0.5944334,-1,0,-1                          , 0,0,0,1,0,0 },
{/*SER.N,   */ 0.792522219,0.507272727,1.008270973,0.624605678,1,1,0                          , 0,0,0,0,0,1 },
{/*SER.CA,  */ 0.704259884,0.578181818,0.831823553,0.651419558,-1,0,-1                        , 0,0,0,1,0,0 },
{/*SER.C,   */ 0.412810297,0.327272727,0.435604569,0.351735016,1,1,0                          , 0,0,0,1,0,0 },
{/*SER.O,   */ 0.541679436,0.533636364,0.51043718,0.441640379,1,1,0                           , 0,0,0,0,1,0 },
{/*SER.CB,  */ 1.267391971,1.083636364,1.659314691,1.186119874,-1,0,-1                        , 0,0,0,1,0,0 },
{/*SER.OG,  */ 2.0452038,1.427272727,2.329657345,1.416403785,1,1,0                            , 0,0,0,0,1,1 },
{/*THR.N,   */ 0.691484185,0.521484375,0.689667374,0.5954323,1,1,0                            , 0,0,0,0,0,1 },
{/*THR.CA,  */ 0.507866991,0.481445313,0.56581741,0.481239804,-1,0,-1                         , 0,0,0,1,0,0 },
{/*THR.C,   */ 0.330413625,0.306640625,0.332271762,0.313213703,1,1,0                          , 0,0,0,1,0,0 },
{/*THR.O,   */ 0.517923763,0.516601563,0.491153574,0.432300163,1,1,0                          , 0,0,0,0,1,0 },
{/*THR.CB,  */ 0.908515815,0.736328125,1.072894551,0.763458401,-1,0,-1                        , 0,0,0,1,0,0 },
{/*THR.OG1, */ 1.797080292,1.455078125,2.029723992,1.414355628,1,1,0                          , 0,0,0,0,1,1 },
{/*THR.CG2, */ 1.00973236,0.625,1.38995046,0.662316476,-1,0,-1                                , 0,0,0,1,0,0 },
{/*TRP.N,   */ 0.199915469,0.358669834,0.221631206,0.198083067,1,1,0                          , 0,0,0,0,0,1 },
{/*TRP.CA,  */ 0.207523246,0.185273159,0.239952719,0.185303514,-1,0,-1                        , 0,0,0,1,0,0 },
{/*TRP.C,   */ 0.196534235,0.154394299,0.229905437,0.134185304,1,1,0                          , 0,0,0,1,0,0 },
{/*TRP.O,   */ 0.33262891,0.204275534,0.363475177,0.392971246,1,1,0                           , 0,0,0,0,1,0 },
{/*TRP.CB,  */ 0.489433643,0.39192399,0.609338061,0.447284345,-1,0,-1                         , 0,0,0,1,0,0 },
{/*TRP.CG,  */ 0.519442096,0.330166271,0.658392435,0.428115016,-1,0,-1                        , 1,0,0,0,0,0 },
{/*TRP.CD1, */ 0.796280642,0.57719715,0.929669031,0.658146965,-1,0,-1                         , 1,0,0,0,0,0 },
{/*TRP.CD2, */ 0.672020287,0.432304038,0.84751773,0.565495208,-1,0,-1                         , 1,0,0,0,0,0 },
{/*TRP.NE1, */ 1.179628064,1.052256532,1.307919622,1.309904153,1,1,0                          , 0,0,0,0,0,1 },
{/*TRP.CE2, */ 0.954353339,0.64608076,1.115248227,0.728434505,-1,0,-1                         , 1,0,0,0,0,0 },
{/*TRP.CE3, */ 0.71386306,0.486935867,0.903664303,0.610223642,-1,0,-1                         , 1,0,0,0,0,0 },
{/*TRP.CZ2, */ 1.1386306,0.712589074,1.313829787,0.769968051,-1,0,-1                          , 1,0,0,0,0,0 },
{/*TRP.CZ3, */ 0.781910397,0.660332542,0.956264775,0.82428115,-1,0,-1                         , 1,0,0,0,0,0 },
{/*TRP.CH2, */ 0.928994083,0.669833729,1.106382979,0.766773163,-1,0,-1                        , 1,0,0,0,0,0 },
{/*TYR.N,   */ 0.204312668,0.334035827,0.253123303,0.162754304,1,1,0                          , 0,0,0,0,0,1 },
{/*TYR.CA,  */ 0.220125786,0.174920969,0.272406301,0.17370892,-1,0,-1                         , 0,0,0,1,0,0 },
{/*TYR.C,   */ 0.168014376,0.197049526,0.197990223,0.17057903,1,1,0                           , 0,0,0,1,0,0 },
{/*TYR.O,   */ 0.273674753,0.210748156,0.308799565,0.297339593,1,1,0                          , 0,0,0,0,1,0 },
{/*TYR.CB,  */ 0.422461815,0.318229715,0.570342205,0.328638498,-1,0,-1                        , 0,0,0,1,0,0 },
{/*TYR.CG,  */ 0.39442947,0.169652266,0.539652363,0.214397496,-1,0,-1                         , 1,0,0,0,0,0 },
{/*TYR.CD1, */ 0.589218329,0.299262381,0.800380228,0.323943662,-1,0,-1                        , 1,0,0,0,0,0 },
{/*TYR.CD2, */ 0.608086253,0.272918862,0.778924498,0.34741784,-1,0,-1                         , 1,0,0,0,0,0 },
{/*TYR.CE1, */ 0.903504043,0.610115911,1.109722977,0.715179969,-1,0,-1                        , 1,0,0,0,0,0 },
{/*TYR.CE2, */ 0.87672956,0.49947313,1.016838675,0.552425665,-1,0,-1                          , 1,0,0,0,0,0 },
{/*TYR.CZ,  */ 0.955256065,0.616438356,1.086094514,0.680751174,-1,0,-1                        , 1,0,0,0,0,0 },
{/*TYR.OH,  */ 2.02821204,1.131717597,2.095871809,1.17370892,1,1,0                            , 0,0,0,0,1,1 },
{/*VAL.N,   */ 0.555642787,0.522207268,0.482459761,0.187116564,1,1,0                          , 0,0,0,0,0,1 },
{/*VAL.CA,  */ 0.312463199,0.238223419,0.341312423,0.208588957,-1,0,-1                        , 0,0,0,1,0,0 },
{/*VAL.C,   */ 0.297350343,0.236877524,0.262484523,0.18404908,1,1,0                           , 0,0,0,1,0,0 },
{/*VAL.O,   */ 0.697154073,0.415881561,0.503920759,0.450920245,1,1,0                          , 0,0,0,0,1,0 },
{/*VAL.CB,  */ 0.386849853,0.259757739,0.485348741,0.242331288,-1,0,-1                        , 0,0,0,1,0,0 },
{/*VAL.CG1, */ 0.947988224,0.63795424,1.221626083,0.858895706,-1,0,-1                         , 0,0,0,1,0,0 },
{/*VAL.CG2, */ 1.067713445,0.63795424,1.494841106,0.763803681,-1,0,-1                         , 0,0,0,1,0,0 }};

//-> feature_selection (according to P2Rank)
const int Atom_FeatDim=9;
int Atom_FeatSelect[Atom_FeatDim]={0,1,4,7,8,9,10,11,12};

//------------- ResiAtom_Mapping -------------//
//ARNDCQEGHILKMFPSTWYVZ
int ResiAtom_Map[Resi_Number]=
{ 0,5,16,24,32,38,47,56,60,70,78,86,95,103,114,121,127,134,148,160};
//A R  N  D  C  Q  E  G  H  I  L  K  M   F   P   S   T   W   Y   V


//=============== residue features =================//
double Resi_Features[Resi_Number][22]={
//selection:      0    0     0     0     0     0          1     1     1     1     1     1    1    1    1     1    1   1    1     1     1    1
//resiName,    pnas1,pnas2,pnas3,pnas4,pnas5, RAx,     hindex,hphob,hphil,aliph,aroma,sulfu,hxyl,base,acid,amide,pos,neg,donor,accep,doacc,ion
{/*A - ALA*/  -0.591,-1.302,-0.733,1.57,-0.146,0.701  ,   1.8, 1,0,1,0,0,0,0,0,0,0,0,0,0,0,0 },
{/*R - ARG*/  1.538,-0.055,1.502,0.44,2.897,0.916     ,  -4.5, 0,1,0,0,0,0,3,0,0,1,0,1,0,0,1 },
{/*N - ASN*/  0.945,0.828,1.299,-0.169,0.933,0.811    ,  -3.5, 0,1,0,0,0,0,0,0,1,0,0,0,0,1,0 },
{/*D - ASP*/  1.05,0.302,-3.656,-0.259,-3.242,1.015   ,  -3.5, 0,1,0,0,0,0,0,1,0,0,1,0,1,0,1 },
{/*C - CYS*/  -1.343,0.465,-0.862,-1.02,-0.255,1.65   ,   2.5, 1,0,0,0,1,0,0,0,0,0,0,0,0,0,1 },
{/*Q - GLN*/  0.931,-0.179,-3.005,-0.503,-1.853,0.669 ,  -3.5, 0,1,0,0,0,0,0,0,1,0,0,0,0,1,0 },
{/*E - GLU*/  1.357,-1.453,1.477,0.113,-0.837,0.956   ,  -3.5, 0,1,0,0,0,0,0,1,0,0,1,0,1,0,1 },
{/*G - GLY*/  -0.384,1.652,1.33,1.045,2.064,0.788     ,  -0.4, 1,0,1,0,0,0,0,0,0,0,0,0,0,0,0 },
{/*H - HIS*/  0.336,-0.417,-1.673,-1.474,-0.078,2.286 ,  -3.2, 0,1,0,0,0,0,1,0,0,1,0,0,0,1,1 },
{/*I - ILE*/  -1.239,-0.547,2.131,0.393,0.816,1.006   ,   4.5, 1,0,1,0,0,0,0,0,0,0,0,0,0,0,0 },
{/*L - LEU*/  -1.019,-0.987,-1.505,1.266,-0.912,1.045 ,   3.8, 1,0,1,0,0,0,0,0,0,0,0,0,0,0,0 },
{/*K - LYS*/  1.831,-0.561,0.533,-0.277,1.648,0.468   ,  -3.9, 0,1,0,0,0,0,2,0,0,1,0,1,0,0,1 },
{/*M - MET*/  -0.663,-1.524,2.219,-1.005,1.212,1.894  ,   1.9, 1,0,0,0,1,0,0,0,0,0,0,0,0,0,0 },
{/*F - PHE*/  -1.006,-0.59,1.891,-0.397,0.412,1.952   ,   2.8, 1,0,0,1,0,0,0,0,0,0,0,0,0,0,0 },
{/*P - PRO*/  0.189,2.081,-1.628,0.421,-1.392,0.212   ,  -1.6, 1,0,1,0,0,0,0,0,0,0,0,0,0,0,0 },
{/*S - SER*/  -0.228,1.399,-4.76,0.67,-2.647,0.883    ,  -0.8, 0,1,0,0,0,1,0,0,0,0,0,0,0,1,0 },
{/*T - THR*/  -0.032,0.326,2.213,0.908,1.313,0.73     ,  -0.7, 0,1,0,0,0,1,0,0,0,0,0,0,0,1,0 },
{/*W - TRP*/  -0.595,0.009,0.672,-2.128,-0.184,3.084  ,  -0.9, 1,0,0,1,0,0,0,0,0,0,0,1,0,0,0 },
{/*Y - TYR*/  0.26,0.83,3.097,-0.838,1.512,1.672      ,  -1.3, 0,1,0,1,0,0,0,0,0,0,0,0,0,1,1 },
{/*V - VAL*/  -1.337,-0.279,-0.544,1.242,-1.262,0.884 ,   4.2, 1,0,1,0,0,0,0,0,0,0,0,0,0,0,0 }};

//-> feature_selection (according to P2Rank)
const int Resi_FeatDim=16;
int Resi_FeatSelect[Resi_FeatDim]={6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};

//-------- char to AminoAcid --------//
//                 A   B   C   D   E   F   G  H    I  J  K  L    M   N   O   P   Q   R    S    T  U   V   W   X   Y   Z
int AA4SUB[26] = { 1, -1,  5,  4,  7, 14,  8, 9,  10,-1, 12,11,  13, 3, -1, 15,  6,  2,  16,  17,-1, 20, 18, -1, 19, -1};



//--------------------- calculate resi/atom-level features ---------------------//
void Calculate_ResiAtom_Features(
	vector <vector <string> > &str,
	vector <vector <vector <double> > > &resiatom_feature)
{
	//-> init
	resiatom_feature.resize(str.size());
	for(int i=0;i<(int)str.size();i++)
	{
		resiatom_feature[i].resize(str[i].size());
		for(int j=0;j<(int)str[i].size();j++)
		{
			resiatom_feature[i][j].resize(Resi_FeatDim+Atom_FeatDim,0);
		}
	}
	//-> get feature
	for(int i=0;i<(int)str.size();i++)
	{
		for(int j=0;j<(int)str[i].size();j++)
		{
			//--| check resi
			char resi_char=str[i][j][0];
			if(resi_char<'A' || resi_char>'Z')continue;
			int resi_int=AA4SUB[resi_char-'A']-1;
			if(resi_int<0)continue;
			//--| check atom
			char atom_char=str[i][j][2];
			if(atom_char<'a' || atom_char=='!')continue;
			int atom_int=atom_char-'a';
			int atom_begi=ResiAtom_Map[resi_int];
			if(atom_begi+atom_int>=Atom_Number)continue;
			int atom_len;
			if(resi_int<Resi_Number-1)atom_len=ResiAtom_Map[resi_int+1]-atom_begi;
			else atom_len=Atom_Number-atom_begi;
			if(atom_int>=atom_len)continue;
			//--| assign resi features
			int count=0;
			for(int k=0;k<Resi_FeatDim;k++)
			{
				resiatom_feature[i][j][count]=Resi_Features[resi_int][Resi_FeatSelect[k]];
				count++;
			}
			//--| assign atom features
			int atom_index=atom_begi+atom_int;
			for(int k=0;k<Atom_FeatDim;k++)
			{
				resiatom_feature[i][j][count]=Atom_Features[atom_index][Atom_FeatSelect[k]];
				count++;
			}
		}
	}
}



//================================== The calculation of Protrusion score =============================//
//-> We use the definition of CX score to calculate protrusion score, according to the following paper:
// CX, an algorithm that identifies protruding atoms in proteins

//----- calculate distance of two points --------//
double Distance2(vector <double> &in1,vector <double> &in2)
{
	double dist2=0;
	for(int i=0;i<3;i++)dist2+=(in1[i]-in2[i])*(in1[i]-in2[i]);
	return dist2;
}
double Distance(vector <double> &in1,vector <double> &in2)
{
	double dist=0;
	for(int i=0;i<3;i++)dist+=(in1[i]-in2[i])*(in1[i]-in2[i]);
	return sqrt(1.0*dist);
}

//---- part 1: get grid from the input XYZ -----//
void Get_Grid_from_XYZ(double grid_interval,
	vector <vector <vector <double> > > &xyz,
	vector <vector <vector <vector <pair<int,int> > > > > &grid)
{
	//-> get min and max coordinate
	double min_x=999999;
	double min_y=999999;
	double min_z=999999;
	double max_x=-999999;
	double max_y=-999999;
	double max_z=-999999;
	for(int i=0;i<(int)xyz.size();i++)
	{
		for(int j=0;j<(int)xyz[i].size();j++)
		{
			//-> get min
			if(xyz[i][j][0]<min_x)min_x=xyz[i][j][0];
			if(xyz[i][j][1]<min_y)min_y=xyz[i][j][1];
			if(xyz[i][j][2]<min_z)min_z=xyz[i][j][2];
			//-> get max
			if(xyz[i][j][0]>max_x)max_x=xyz[i][j][0];
			if(xyz[i][j][1]>max_y)max_y=xyz[i][j][1];
			if(xyz[i][j][2]>max_z)max_z=xyz[i][j][2];
		}
	}
	int main_x=(int)( 1.0*(max_x-min_x)/grid_interval)+1;
	int main_y=(int)( 1.0*(max_y-min_y)/grid_interval)+1;
	int main_z=(int)( 1.0*(max_z-min_z)/grid_interval)+1;
	//-> create grid
	grid.clear();
	grid.resize(main_x);
	for(int i=0;i<main_x;i++)
	{
		grid[i].resize(main_y);
		for(int j=0;j<main_y;j++)
		{
			grid[i][j].resize(main_z);
		}
	}
	//-> fill-in points
	for(int i=0;i<(int)xyz.size();i++)
	{
		for(int j=0;j<(int)xyz[i].size();j++)
		{
			int index_x=(int)(1.0*(xyz[i][j][0]-min_x)/grid_interval);
			int index_y=(int)(1.0*(xyz[i][j][1]-min_y)/grid_interval);
			int index_z=(int)(1.0*(xyz[i][j][2]-min_z)/grid_interval);
			grid[index_x][index_y][index_z].push_back(pair<int,int>(i,j));
		}
	}
}

//---------- calculate protrusion score -------------//
//-> extract candidate points
void Extract_Candidate_Points(
	int x,int y,int z,int grid_range,
	vector <vector <vector <vector <pair<int,int> > > > > &grid,
	vector <pair<int,int> > &out)
{
	out.clear();
	for(int i=x-grid_range;i<=x+grid_range;i++)
	{
		if(i<0)continue;
		if(i>=grid.size())break;
		for(int j=y-grid_range;j<=y+grid_range;j++)
		{
			if(j<0)continue;
			if(j>=grid[i].size())break;
			for(int k=z-grid_range;k<=z+grid_range;k++)
			{
				if(k<0)continue;
				if(k>=grid[i][j].size())break;
				//-- collect --//
				for(int l=0;l<(int)grid[i][j][k].size();l++)out.push_back(grid[i][j][k][l]);
				
			}
		}
	}
}
//-> calculate protrusion inner part
double Calculate_Protrusion_inner(
	vector <double> &in, double radius,
	vector <pair<int,int> > &candidate,
	vector <vector <vector <double> > > &xyz)
{
	//init
	double rad2=radius*radius;
	double Vsphere=4.0/3.0*3.1415926*radius*radius*radius;
	double Vatom=20.1;
	//calculate
	int natom=0;
	for(int k=0;k<(int)candidate.size();k++)
	{
		int first=candidate[k].first;
		int second=candidate[k].second;
		double dist2=Distance2(in,xyz[first][second]);
		if(dist2<rad2)natom++;
	}
	//protrusion score
	double Vint=1.0*natom*Vatom;
	double Vext=Vsphere-Vint;
	if(Vext<0)Vext=0;
	return 1.0*Vext/Vint;
}

//------- calculate protrusion main part -------//
void Calculate_Protrusion(
	double radius,double grid_interval,
	vector <vector <vector <double> > > &xyz,
	vector <vector <double> > &protrusion_score)
{
	//-> init
	protrusion_score.resize(xyz.size());
	for(int i=0;i<(int)xyz.size();i++)protrusion_score[i].resize(xyz[i].size());
	//-> calculate grid	
	int grid_range=1.0*radius/grid_interval+1;
	vector <vector <vector <vector <pair<int,int> > > > > grid;
	Get_Grid_from_XYZ(grid_interval,xyz,grid);
	//-> traverse all cells in the grid
	for(int i=0;i<(int)grid.size();i++)
	{
		for(int j=0;j<(int)grid[i].size();j++)
		{
			for(int k=0;k<(int)grid[i][j].size();k++)
			{
				//--| get point size
				int size=(int)grid[i][j][k].size();
				if(size<=0)continue;
				//--| generate candidate for each cell
				vector <pair<int,int> > candidate;
				Extract_Candidate_Points(i,j,k,grid_range,grid,candidate);
				//--| for each point
				for(int l=0;l<size;l++)
				{
					int first=grid[i][j][k][l].first;
					int second=grid[i][j][k][l].second;
					vector <double> point=xyz[first][second];
					double score=Calculate_Protrusion_inner(point,radius,candidate,xyz);
					protrusion_score[first][second]=score;
				}
			}
		}
	}
}


//========= XYZ format for point-cloud ==========//
/*
    6 DOg  -22.389   29.406   56.176 0  58
    6 DOh  -23.720   31.160   56.319 0  58
    7 KNa  -22.679   27.076   60.272 0  55
    7 KCb  -21.836   26.578   61.340 0  55
    7 KCc  -20.356   26.635   60.950 0  55
    7 KOd  -19.746   25.631   60.582 0  55
    7 KCe  -22.307   25.189   61.783 0  55
    7 KCf  -23.751   25.237   62.298 0  55
    7 KCg  -24.043   24.207   63.371 0  55
    7 KCh  -25.357   24.538   64.095 0  55
    7 KNi  -25.610   23.600   65.237 0  55
    8 QNa  -19.785   27.829   61.073 0  24
    8 QCb  -18.384   28.110   60.732 0  24
    8 QCc  -17.293   27.383   61.502 0  24
    8 QOd  -16.145   27.388   61.063 0  24
...
*/
//----- check insert_code ------//
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

//--------- load XYZ ----------//
void Load_XYZ_Label(string &fn,vector <vector <vector <double> > > &xyz,
	vector <vector <string> > &str, vector <vector <string> > &lab,
	vector <vector <string> > &remain, vector <string> &resi)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(fn.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",fn.c_str());
		exit(-1);
	}
	xyz.clear();
	str.clear();
	lab.clear();
	remain.clear();
	resi.clear();
	vector <double> point(3);
	vector <vector <double> > xyz_tmp;
	vector <string> str_tmp;
	vector <string> lab_tmp;
	vector <string> remain_tmp;
	string prev="";
	string str_rec;
	int first=1;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		www>>temp;
		//check ins_code
		if(Check_Ins(temp)!=1)temp.push_back(' ');
		//check chain_code
		if(Check_Chain(temp)==0)temp=" "+temp;
		//record first
		if(first==1)
		{
			first=0;
			prev=temp;
		}
		if(temp!=prev)
		{
			xyz.push_back(xyz_tmp);
			str.push_back(str_tmp);
			lab.push_back(lab_tmp);
			remain.push_back(remain_tmp);
			resi.push_back(prev);
			xyz_tmp.clear();
			str_tmp.clear();
			lab_tmp.clear();
			remain_tmp.clear();
			prev=temp;
		}
		www>>temp;
		str_tmp.push_back(temp);
		www>>point[0]>>point[1]>>point[2];
		xyz_tmp.push_back(point);
		//remain
		string remain_rec="";
		int count=0;
		string lab_value="0";
		for(;;)
		{
			if(! (www>>temp) )break;
			if(count>0)remain_rec=remain_rec+temp+" ";
			else lab_value=temp;
			count++;
		}
		lab_tmp.push_back(lab_value);
		remain_tmp.push_back(remain_rec);
	}
	//termi
	if(first==0)
	{
		xyz.push_back(xyz_tmp);
		str.push_back(str_tmp);
		lab.push_back(lab_tmp);
		remain.push_back(remain_tmp);
		resi.push_back(prev);
	}
}



//==================== XYZ Generate Features ======================//
void XYZ_Gen_Feature(string &in, string &out,
	double radius,double grid_interval)
{
	//----- load XYZ -----//
	vector <vector <vector <double> > > xyz;
	vector <vector <string> > str;
	vector <vector <string> > lab;
	vector <vector <string> > remain;
	vector <string> resi;
	Load_XYZ_Label(in,xyz,str,lab,remain,resi);
	//----- generate resi/atom features ----//
	vector <vector <vector <double> > > resiatom_feature;
	Calculate_ResiAtom_Features(str,resiatom_feature);
	//----- generate protrusion features ---//
	vector <vector <double> > protrusion_score;
	Calculate_Protrusion(radius,grid_interval,xyz,protrusion_score);
	//----- output XYZ ------//
	FILE *fp=fopen(out.c_str(),"wb");
	for(int i=0;i<(int)xyz.size();i++)
	{
		for(int j=0;j<(int)xyz[i].size();j++)
		{
			//transfer resi/atom feature to string
			stringstream oss;
			for(int k=0;k<(int)resiatom_feature[i][j].size();k++)
			{
				int wsiii=(int)resiatom_feature[i][j][k];
				if(wsiii!=resiatom_feature[i][j][k])oss << resiatom_feature[i][j][k] << " ";
				else oss << wsiii << " ";
			}
			string wsbuf=oss.str();
			//output to file
			fprintf(fp,"%7s %3s %8.3f %8.3f %8.3f %s %s %5.2f %s\n",
				resi[i].c_str(),str[i][j].c_str(),
				xyz[i][j][0],xyz[i][j][1],xyz[i][j][2],
				lab[i][j].c_str(),remain[i][j].c_str(),
				protrusion_score[i][j],wsbuf.c_str());
		}
	}
	fclose(fp);
}


//-------- main ---------//
int main(int argc,char **argv)
{
	//------ XYZ_Gen_Feature -------//
	{
		if(argc<5)
		{
			fprintf(stderr,"XYZ_Gen_Feature <xyz_in> <xyz_out> <radius> <grid_interval> \n");
			fprintf(stderr,"[note]: set radius to 10; set grid_interval to 1 or 2;\n");
			exit(-1);
		}
		string xyz_in=argv[1];
		string xyz_out=argv[2];
		double radius=atof(argv[3]);
		double grid_interval=atof(argv[4]);
		//proc
		XYZ_Gen_Feature(xyz_in,xyz_out,radius,grid_interval);
		//exit
		exit(0);
	}
}



