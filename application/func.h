#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <queue>
#include <cstdlib>
#include <time.h>
#include <float.h>
#include <sstream>

using namespace std;

int initialize(void);
int reInitialize(void);
int initializeEvent(void);
int readInput(int, char**, int);
int checkFirstPrim(int);
int checkFirstPrimRerouting(int);
int retuneDownRerouting();
int asign(int, int);
int asignRerouting(int,int);
int asignBp(int, int);
int asignBpRerouting(int, int);
int printSpec(void);
int genDemands(int);
void delFromList(int, int);
void addToList(int, int, int);
int checkExactPrim(int);
int checkExactBack(int);
int checkExactPrimRerouting(int);
int checkExactBackRerouting(int);
int searchPrimRoute(int, int);
int searchBackRoute(int, int);
void deleteLPRerouting(int, int);
int removeLP1_1(int);
int firstFit1_1(int);
int checkFirstBack(int);
int checkFirstBackRerouting(int);
int retuneBp(int);
int readResultPy();
int readResultReroutingPy();
void statDefragPy(int);
void statDefragReroutingPy(int);
int statAlgoRerouting(void);
int writeOutput(int);
int writeOutputPy(int);
int writeOutputReroutingPy(int);
void addToList2(int, int, int);
void delFromList2(int, int, int);
double sum(double, int);
double ave(double, int);
double var(double, int);
double standard(double, int);
double finTime();

int lp_size[REQ_NUM], source[REQ_NUM], dest[REQ_NUM], spec_ind[REQ_NUM];
bool spec[CAPASITY][LINK_NUM];
bool path_prim[NODE_NUM][NODE_NUM][LINK_NUM];
bool path_back[NODE_NUM][NODE_NUM][LINK_NUM];
bool linked_path[NODE_NUM][NODE_NUM][NODE_NUM][NODE_NUM];
bool linked_bp[NODE_NUM][NODE_NUM][NODE_NUM][NODE_NUM];
bool linked_crosspath[NODE_NUM][NODE_NUM][NODE_NUM][NODE_NUM];
int blocked;
int isactive[REQ_NUM];
double t_req_event[REQ_NUM], t_hold_event[REQ_NUM], t_exp_event[REQ_NUM];
int hops[NODE_NUM][NODE_NUM], bhops[NODE_NUM][NODE_NUM];
int link[NODE_NUM][NODE_NUM];
int last_lp;
int midlp = REQ_NUM;
unsigned seed1 = 1448601515;
unsigned seed2 = 125;
struct lpNode{
	int x;
	int y;
	int z;
	lpNode *next;
};
lpNode *activeList = NULL;
lpNode *backupList = NULL;
lpNode *mixtList = NULL;
lpNode *realList = NULL;
lpNode *tempList = NULL;
struct Node {
	vector<int> edges_to;
 	vector<int> edges_cost;
	bool done;
	int cost;
	int nodeNum;
	int from;
	bool operator> (const Node &node) const {
   	return (cost > node.cost);
	}
};
bool path_prim_rr[LINK_NUM][REQ_NUM], path_back_rr[LINK_NUM][REQ_NUM];
int part[NODE_NUM][NODE_NUM];
int algoCall;
int t_temp =0;
int last_ret= 0;
double time_slot_now;
int bp_ind[REQ_NUM];
int togOp;
int realOp;
int rerouteOp;
bool lpState[REQ_NUM];
bool bpState[REQ_NUM];
int last_blocked = 0;
struct Event
{
	double time_slot_now;
	int lpNum;
	int type;	
	bool operator> (const Event &event) const {
   	return (time_slot_now > event.time_slot_now);
	}
};
priority_queue<Event, vector<Event>, greater<Event> > eventQueue;
priority_queue<Event, vector<Event>, greater<Event> > deleteQueue;
Event startEvent[REQ_NUM];
Event endEvent[REQ_NUM];
Event nowEvent;
Event nextEvent;
int defragCount;
vector<Event> defragEvent;

int readInput(int argc, char* argv[0], int load)
{
	ifstream fin;
	// register network
	stringstream ss;
	string s;
	ss << "./../network/net_" << argv[1] << ".txt";
	s = ss.str();
	fin.open (s);
	if (!fin){
		cout <<"[error] cannot open network model" << endl;
		return 1;
	}
	fin.ignore(INT_MAX,'=');
	for(int k = 0; k < LINK_NUM; k++)
	{
		int a, b;
		fin >> a >> b ;
		link[a][b]=k;
	}
	fin.close();

	// count primary path number
	ss.str("");
	ss << "./../network/path_" << argv[1] << "_prim.txt";
	s = ss.str();
	fin.open (s);
	if (!fin)
	{
		cout <<"[error] cannot open primary routing file" << endl;
		return 1;
	}
	char tmp;
	fin.ignore(INT_MAX,'=');
	fin >> tmp;
	int path_num = 0;
	while(tmp !=';')
	{
		path_num++;
		fin.ignore(INT_MAX,'\n');
		fin >> tmp;
	}
	fin.close();

	// registar primary path
	ss.str("");
	ss << "./../network/path_" << argv[1] << "_prim.txt";
	s = ss.str();
	fin.open (s);
	fin.ignore(INT_MAX,'=');
	for (int k = 0; k < path_num; k++)
	{
		int i, j, a, b, p;
		fin >> i >> j >> a >> b;
		fin.ignore(INT_MAX,'\n');
		p = link[a][b];
		path_prim[i][j][p] = 1;
		++hops[i][j];
	}
	fin.close();

	// count backup path number
	ss.str("");
	ss << "./../network/path_" << argv[1] << "_back.txt";
	s = ss.str();
	fin.open (s);
	if (!fin)
	{
		cout <<"[error] cannot open backup routing file" << endl;
		return 1;
	}
	fin.ignore(INT_MAX,'=');
	fin >> tmp;
	path_num = 0;
	while(tmp !=';'){
		path_num++;
		fin.ignore(INT_MAX,'\n');
		fin >> tmp;
	}
	fin.close();

	// registar backup path
	ss.str("");
	ss << "./../network/path_" << argv[1] << "_back.txt";
	s = ss.str();
	fin.open (s);
	fin.ignore(INT_MAX,'=');
	for (int k=0; k < path_num; k++)
	{
		int i, j, a, b, p;
		fin >> i >> j >> a >> b;
		fin.ignore(INT_MAX,'\n');
		p = link[a][b];
		path_back[i][j][p] = 1;
		++bhops[i][j];
	}
	fin.close();

	for (int i = 0; i < NODE_NUM; i++)
	{
		for (int j = 0; j < NODE_NUM; j++)
		{
			for (int k = 0; k < NODE_NUM; k++)
			{
				for (int l = 0; l < NODE_NUM; l++)
				{
					if (k == i && l == j){
						linked_path[i][j][k][l] = 0;
					} else {
						for(int p = 0; p < LINK_NUM; p++){
							if(path_prim[i][j][p] && path_prim[k][l][p]) linked_path[i][j][k][l] = 1;
							if(path_back[i][j][p] && path_back[k][l][p]) linked_bp[i][j][k][l] = 1;
							if(path_prim[i][j][p] && path_back[k][l][p]) linked_crosspath[i][j][k][l] = 1;
						}
					}
				}
			}
		}
	}
	return 0;
}

int writeOutput(int load)
{
	int i,j,ind=0;
	int lp,n,f0;
	int s,d;
	int c = MAX_STEP;
	int m=0;

	ofstream ofs2;
	ofs2.open("smpe_output.dat");
	if(!ofs2){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}

	lpNode *cur = activeList;
	while ( cur != NULL ) {						// Checking all active LPs
		m++;
		lp = cur->x;
		cur = cur->next;
	}

	ofs2 << endl;
	ofs2 << "param REQ_NUM := " << m <<";" << endl;
	ofs2 << "param B := " << CAPASITY <<";" << endl;
	ofs2 << "param LINK_NUM := " << LINK_NUM <<";" << endl;
	ofs2 << "param C := " << c <<";" << endl;
	ofs2 << endl;
	ofs2 << "param : CAPASITY K NODE_NUM F0 := " << endl;

	ofstream ofs3;
	ofs3.open("running.txt");
	if(!ofs3){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}
	ofs3 << "Load load := " << load <<";" << endl;
	ofs3 << "param REQ_NUM := " << m <<";" << endl;
	ofs3 << "Last LP := " << lp <<";" << endl;
	ofs3.close();
	// cout << "Load load := " << load <<";" << endl;
	// cout << "param REQ_NUM := " << m <<";" << endl;
	// cout << "Last LP := " << lp <<";" << endl;

	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		n = lp_size[lp];
		f0 = spec_ind[lp];
		ofs2 << 2*ind <<" "<< ind <<" 1 "<< n <<" "<< f0 << endl;
//	if(ind<3 || ind> 27) cout << "LP " << lp << " at " << f0 << endl;
		f0 = bp_ind[lp];
		ofs2 << 2*ind+1 <<" "<< ind <<" 0 "<< n <<" "<< f0 << endl;
//	if(ind<3 || ind> 27) cout << "LP " << lp << " at " << f0 << endl;
		ind++;
	}
	ofs2 << ";" << endl << endl;
	cout << endl;

	ofs2 << "set Q := " << endl;
	ind =0;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		s = source[lp];
		d= dest[lp];

		for(j=0;j<LINK_NUM;j++){
			if(path_prim[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
		for(j=0;j<LINK_NUM;j++){
			if(path_back[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
	}
	ofs2 << ";" << endl << endl;

	ofs2.close();
	return 0;
}

int writeOutputPy(int load)
{
	int i,j,ind=0;
	int lp,n,f0,k;
	int s,d;
	int c = MAX_STEP;
	int m=0;

	ofstream ofs2;
	ofs2.open("./../result/ssr_lno_output.txt");
	if(!ofs2){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}

	lpNode *cur = activeList;
	while ( cur != NULL ) {						// Checking all active LPs
		m++;   //activeなリクエストの数を数える
		lp = cur->x;
		cur = cur->next;
	}

	ofs2 << "param CAPASITY   := " << m   << endl;
	ofs2 << "param absP:= " << m*2 << endl;
	ofs2 << "param absF:= " << CAPASITY   << endl;
	ofs2 << "param absE:= " << LINK_NUM   << endl;
	ofs2 << "param T   := " << c   << endl;
	ofs2 << endl;

	ofstream ofs3;
	ofs3.open("./../result/ssr_lno_running.txt");
	if(!ofs3){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}
	ofs3 << "Load load := " << load <<";" << endl;
	ofs3 << "param REQ_NUM := " << m <<";" << endl;
	ofs3 << "Last LP := " << lp <<";" << endl;
	ofs3.close();
	// cout << "Load load := " << load <<";" << endl;
	// cout << "param REQ_NUM := " << m <<";" << endl;
	// cout << "Last LP := " << lp <<";" << endl;


	ofs2 << "param : s_p k_p_init n_p f_p_init := " << endl;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		n = lp_size[lp];
		k = lpState[lp];
		if (k)
		{
			f0 = spec_ind[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			f0 = bp_ind[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			ind++;
		}else{
			f0 = spec_ind[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			f0 = bp_ind[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			ind++;
		}
	}
	ofs2 << ";" << endl << endl;
	// cout << endl;

	ofs2 << "set P_e := " << endl;
	ind =0;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		s = source[lp];
		d= dest[lp];

		for(j=0;j<LINK_NUM;j++){
			if(path_prim[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
		for(j=0;j<LINK_NUM;j++){
			if(path_back[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
	}
	ofs2 << ";" << endl << endl;

	// printSpec();

	ofs2.close();
	return 0;
}

int writeOutputReroutingPy(int load)
{
	int i,j,ind=0;
	int lp,n,f0,k;
	int s,d;
	int c = MAX_STEP;
	int m=0;

	ofstream ofs2;
	ofs2.open("./../result/ssrr_lno_output.txt");
	if(!ofs2){
		cout<< "Cannot open ssrr_lno_output.txt file"<<endl;
		return 1;
	}

	lpNode *cur = activeList;
	while ( cur != NULL ) {						// Checking all active LPs
		m++;   //activeなリクエストの数を数える
		lp = cur->x;
		cur = cur->next;
	}

	ofs2 << "param CAPASITY   := " << m   << endl;
	ofs2 << "param absP:= " << m*2 << endl;
	ofs2 << "param absF:= " << CAPASITY   << endl;
	ofs2 << "param absE:= " << LINK_NUM   << endl;
	ofs2 << "param absV:= " << NODE_NUM   << endl;
	ofs2 << "param T   := " << c   << endl;
	ofs2 << endl;

	ofstream ofs3;
	ofs3.open("./../result/ssr_lno_running.txt");
	if(!ofs3){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}
	ofs3 << "Load load := " << load <<";" << endl;
	ofs3 << "param REQ_NUM := " << m <<";" << endl;
	ofs3 << "Last LP := " << lp <<";" << endl;
	ofs3.close();
	// cout << "Load load := " << load <<";" << endl;
	// cout << "param REQ_NUM := " << m <<";" << endl;
	// cout << "Last LP := " << lp <<";" << endl;


	ofs2 << "param : s_p k_p_init n_p f_p_init := " << endl;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		n = lp_size[lp];
		k = lpState[lp];
		if (k)
		{
			f0 = spec_ind[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			f0 = bp_ind[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			ind++;
		}else{
			f0 = spec_ind[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			f0 = bp_ind[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			ind++;
		}
	}
	ofs2 << ";" << endl << endl;
	// cout << endl;

	ofs2 << "E := " << endl;
	for ( i = 0; i < NODE_NUM; i++){
		for ( j = 0; j < NODE_NUM; j++){
			if (link[i][j] < LINK_NUM){
				ofs2 << i << " " << j << endl;
			}
		}
	}
	ofs2 << ";" << endl << endl;

	ofs2 << "r_p_i_j := " << endl;
	// cout << "r_p_i_j := " << endl;
	ind =0;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		for ( i = 0; i < NODE_NUM; i++) {
			for ( j = 0; j < NODE_NUM; j++ ){
				if(link[i][j] < LINK_NUM){
					if (path_prim_rr[link[i][j]][lp]){
						ofs2 << ind << " " << i << " " << j << endl;
						// cout << ind << " " << i << " " << j << endl;
					}
				}
			}
		}

		ind++;//バックアップパスが連番になっている
		for ( i = 0; i < NODE_NUM; i++) {
			for ( j = 0; j < NODE_NUM; j++ ){
				if(link[i][j] < LINK_NUM){
					if (path_back_rr[link[i][j]][lp]){
						ofs2 << ind << " " << i << " " << j << endl;
						// cout << ind << " " << i << " " << j << endl;
					}
				}
			}
		}
		ind++;
	}
	ofs2 << ";" << endl << endl;

	ofs2 << "o_p d_p := " << endl;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		s = source[lp];
		d= dest[lp];
		ofs2 << s << " " << d << endl;
	}
	ofs2 << ";" << endl << endl;

	// printSpec();

	ofs2.close();
	return 0;
}


int retuneBp(int load)
{
	if(algoCall == 0)
	{
		 return 0;
	}

	if(algoCall==1){
		statAlgoRerouting();
		retuneDownRerouting();
		return 0;
	}

	if(algoCall==2){
		statAlgoRerouting();
		retuneDownRerouting();
		return 0;
	}

	if(algoCall==3){
		statDefragPy(load);
		return 0;
	}

	if(algoCall==4){
		statDefragReroutingPy(load);
		return 0;
	}
	return 0;
}

double sum(double data[], int n)
{
    int i;                  // 変数の宣言
    double total = 0.0;
   
    for(i=0; i<n; i++){
        total += data[i];   // 合計を計算
    }
       
    return total;           // 合計を返す
}
 
// 平均の計算
double ave(double data[], int n)
{
    double total = sum(data, n);    // 合計
    return total/n;                     // 平均=合計/データ数
}

// 分散の計算
double var(double data[], int n) {
    int i;
    double a = ave(data, n);    // 平均値
    double v = 0.0;                             // 分散
    // 分散を計算
    for (i=0; i<n; i++)
        v += (data[i] - a) * (data[i] - a);
    return v/n;                                 // 分散の平均
}

// 標準偏差の計算
double standard(double data[], int n) {
    return sqrt(var(data, n));                  // 標準偏差=分散の平方根
}

void statDefragPy(int load)
{
	writeOutputPy(load);
	system("python ssr_lno.py");
	readResultPy();
}

void statDefragReroutingPy(int load)
{
	writeOutputReroutingPy(load);
	system("python ssrr_lno.py");
	readResultReroutingPy();
}

int readResultPy()
{
	int i,j;
	int a,b, lp;
	ifstream fin, fin1;

	for(i=0;i<CAPASITY;i++){
		for(j=0;j<LINK_NUM;j++)  spec[i][j]= 0;
	}

	fin.open ("./../result/ssr_lno_result.txt");
		if (!fin){
			cout <<"Cannot read ssr_lno_result.txt file" << endl;
			return 1;
		}

		fin.ignore(INT_MAX,'=');
		fin >> a ;
		togOp += a;

		fin.ignore(INT_MAX,'=');
		fin >> b ;
		realOp += b;

		fin.ignore(INT_MAX,'=');

		lpNode *cur = activeList;
		while ( cur != NULL ) {			// Same as 		for(k=0; k<m; k++){ with m nb of siganls
			lp = cur->x;
			cur = cur->next;

			fin >> a >> b ;
			asign(lp, b);
			fin.ignore(INT_MAX,'\n');
			fin >> a >> b ;
			asignBp(lp, b);
			fin.ignore(INT_MAX,'\n');
		}
	fin.close();
	return 0;
}

int readResultReroutingPy()
{
	int i,j;
	int a,b,c,d,lp;
	ifstream fin, fin1;

	for(i=0;i<CAPASITY;i++){
		for(j=0;j<LINK_NUM;j++)  spec[i][j]= 0;
	}

	fin.open ("./../result/ssrr_lno_result.txt");
		if (!fin){
			cout <<"Cannot read ssrr_lno_result.txt file" << endl;
			return 1;
		}

		fin.ignore(INT_MAX,'=');
		fin >> a ;
		togOp += a;

		fin.ignore(INT_MAX,'=');
		fin >> b ;
		realOp += b;

		lpNode *cur = activeList;
		while ( cur != NULL ) {			// Same as 		for(k=0; k<m; k++){ with m nb of siganls
			lp = cur->x;
			cur = cur->next;

			//initialize
			for(i=0;i<LINK_NUM;i++){
				path_prim_rr[i][lp] = 0;
				path_back_rr[i][lp]   = 0;
			}

			fin.ignore(INT_MAX,'=');
			fin >> a >> b >> c >> d;
			while(!d){//同じパスなら
				path_prim_rr[link[b][c]][lp] = 1;
				fin.ignore(INT_MAX,'\n');
				// cout << "path_prim_rr[link[" << b << "][" << c << "]][" << lp << "] = " << path_prim_rr[link[b][c]][lp] << endl;
				fin >> a >> b >> c >> d;
			}
			fin.ignore(INT_MAX,'=');
			fin >> a >> b >> c >> d;
			while(!d){//同じパスなら
				path_back_rr[link[b][c]][lp] = 1;
				fin.ignore(INT_MAX,'\n');
				// cout << "path_back_rr[link[" << b << "][" << c << "]][" << lp << "] = " << path_back_rr[link[b][c]][lp] << endl;
				fin >> a >> b >> c >> d;
			}
		}

		fin.ignore(INT_MAX,'=');

		cur = activeList;
		while ( cur != NULL ) {			// Same as 		for(k=0; k<m; k++){ with m nb of siganls
			lp = cur->x;
			cur = cur->next;

			fin >> a >> b ;
			// printSpec();
			asignRerouting(lp, b);
			fin.ignore(INT_MAX,'\n');
			fin >> a >> b ;
			asignBpRerouting(lp, b);
			fin.ignore(INT_MAX,'\n');
		}
	fin.close();
	return 0;
}

int statAlgoRerouting()
{
	int a, b, lp, lp2;
	double ret_time = 0;
	int realMov = 1;
	int s1, d1, s2, d2;
	int st;
	int i, j;

	lpNode *cur;
	lpNode *cur1;
	lpNode *cur2;

	tempList = NULL;
	cur = mixtList;
	while ( cur != NULL ) {
		addToList2(2, cur->x, cur->z);
		cur= cur->next;
	}
	cur = tempList;
	double fin_time;
	fin_time = finTime();
	if (fin_time < 0)
	{
		return 0;
	}
	while(realMov){
		realMov = 0;
		tempList = NULL;
		cur = mixtList;
		while ( cur != NULL ) {
			addToList2(2, cur->x, cur->z);
			cur= cur->next;
		}
		cur = tempList;
		realList = NULL;
		cur1 = NULL;
		while ( cur != NULL ){
			realList = NULL;
			cur1 = cur;
			cur = cur->next;
			while ( cur1 != NULL ) {
				lp = cur1->x;
				st = cur1->z;
				cur1 = cur1->next;
				cur2 = realList;
				if(cur2 == NULL){
					addToList2(1, lp, st);
					if(cur && cur->x == lp && cur->z == st) cur = cur->next;
					delFromList2(4, lp, st);
				}else{
					int conf = 0;
					while(cur2 != NULL){
						if(lpState[cur2->x]){
							if(lpState[lp]){
								for(i=0;i<LINK_NUM;i++){
									if(path_prim_rr[i][cur2->x] && path_prim_rr[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<LINK_NUM;i++){
									if(path_prim_rr[i][cur2->x] && path_back_rr[i][lp])conf =1;
								}

							}
						}else{
							if(lpState[lp]){
								for(i=0;i<LINK_NUM;i++){
									if(path_back_rr[i][cur2->x] && path_prim_rr[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<LINK_NUM;i++){
									if(path_back_rr[i][cur2->x] && path_back_rr[i][lp]) conf =1;
								}
							}
						}
						if(conf) break;
						cur2 = cur2->next;
					}
					if(conf==0){
						addToList2(1, lp, st);
						if(cur && cur->x == lp && cur->z == st) cur = cur->next;
						delFromList2(4, lp, st);
					}
				}
			}
			cur2 = realList;
			int realcheck = realOp;
			while ( cur2 != NULL ){
				lp = cur2->x;
				b =  lpState[lp];
				cur2 = cur2->next;
				if(b){
					int path_rr_prev[LINK_NUM];
					for (int i = 0; i < LINK_NUM; ++i)
					{
						path_rr_prev[i] = 0;
					}
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							path_rr_prev[link[i][j]] = path_prim_rr[link[i][j]][lp];
						}
					}
					deleteLPRerouting(lp, 1);
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkExactPrimRerouting(lp);
					} else {
						a = checkExactPrim(lp);
					}
					if (a == CAPASITY || a > spec_ind[lp])
					{
						if (algoCall == 2 || algoCall == 4)
						{
							a = checkFirstPrimRerouting(lp);
						} else {
							a = checkFirstPrim(lp);
						}
					}
					bool isSame = true;
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							if (path_rr_prev[link[i][j]] != path_prim_rr[link[i][j]][lp])
							{
								isSame = false;
							}
						}
					}
					if (isSame == false)
					{
						rerouteOp++;
					}
					if(a != spec_ind[lp]){
						realOp++;
						realMov++;
						if(lpState[lp]){
							bpState[lp] = lpState[lp];
							lpState[lp] = !bpState[lp];
							togOp++;
						}
					}
					asignRerouting(lp, a);
				}
				if(!b){
					int bp_rr_prev[LINK_NUM];
					for (int i = 0; i < LINK_NUM; ++i)
					{
						bp_rr_prev[i] = 0;
					}
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							bp_rr_prev[link[i][j]] = path_back_rr[link[i][j]][lp];
						}
					}
					deleteLPRerouting(lp, 2);
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkExactBackRerouting(lp);
					} else {
						a = checkExactBack(lp);
					}
					if (a == CAPASITY || a > bp_ind[lp])
					{
						if (algoCall == 2 || algoCall == 4)
						{
							a = checkFirstBackRerouting(lp);
						} else {
							a = checkFirstBack(lp);
						}
					}
					bool isSame = true;
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							if (bp_rr_prev[link[i][j]] != path_back_rr[link[i][j]][lp])
							{
								isSame = false;
							}
						}
					}
					if (isSame == false)
					{
						rerouteOp++;
					}
					if(a != bp_ind[lp]){
						realOp++;
						realMov++;
						if(bpState[lp]){
							bpState[lp] = lpState[lp];
							lpState[lp] = !bpState[lp];
							togOp++;
						}
					}
					asignBpRerouting(lp, a);
				}
			}
			cur2 = realList;
			while ( cur2 != NULL ){
				lp = cur2->x;
				cur2 = cur2->next;
				delFromList(3, lp);
			}
			if(realcheck != realOp){
				ret_time += DEFRAG_TIME;
				if((time_slot_now+ret_time >= fin_time || ret_time >= DEFRAG_TOTAL_TIME_MAX) || eventQueue.empty()){
					return 0;	
				}
			}
		}
	}
	return 0;
}

double finTime()
{
	while(eventQueue.top().type != 0){
		if (eventQueue.empty())
		{
			return -1;
		}
		if (eventQueue.top().type == 1)
		{
			deleteQueue.push(eventQueue.top());
		}
		eventQueue.pop();
	}
	return eventQueue.top().time_slot_now;
}

int retuneDownRerouting()
{
	int s, d;
	int index1;
	int i,j,p;
	int k=0;
	int a;

	int ret_time = 0 ;
	int mov_time = 0 ;

	index1 = 1;

	while(index1 < CAPASITY-1 && k++ < INF){			//Sweep the spectrum
		lpNode *cur = activeList; 
		mov_time = 0 ;
		while ( cur != NULL ) {						// Checking all active LPs
			i = cur->x;
			cur = cur->next;
			if(isactive[i] == 1 && spec_ind[i] == index1){
				int path_rr_prev[LINK_NUM]; //to compare the previous root and new route
				for (int b = 0; b < LINK_NUM; ++b)
				{
					path_rr_prev[b] = 0;
				}
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						path_rr_prev[link[b][c]] = path_prim_rr[link[b][c]][i];
					}
				}
				deleteLPRerouting(i, 1);
				if (algoCall == 2 || algoCall == 4)
				{
					a = checkExactPrimRerouting(i);
				} else {
					a = checkExactPrim(i);
				}
				if(a==CAPASITY || a > spec_ind[i]){
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkFirstPrimRerouting(i);
					} else {
						a = checkFirstPrim(i);
					}
				}
				bool isSame = true;
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						if (path_rr_prev[link[b][c]] != path_prim_rr[link[b][c]][i])
						{
							isSame = false;
						}
					}
				}
				if (isSame == false)
				{
					rerouteOp++;
				}
				if(a != spec_ind[i]){
					realOp++;
					if(lpState[i]){
						bpState[i] = lpState[i];
						lpState[i] = !bpState[i];
						togOp++;
					}
					mov_time++;
				}
				asignRerouting(i, a);
			}

			if(isactive[i] == 1 && bp_ind[i] == index1){
				int bp_rr_prev[LINK_NUM]; //to compare the previous root and new route
				for (int b = 0; b < LINK_NUM; ++b)
				{
					bp_rr_prev[b] = 0;
				}
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						bp_rr_prev[link[b][c]] = path_back_rr[link[b][c]][i];
					}
				}
				deleteLPRerouting(i, 2);
				if (algoCall == 2 || algoCall == 4)
				{
					a = checkExactBackRerouting(i);
				} else {
					a = checkExactBack(i);
				}
				if(a==CAPASITY || a > bp_ind[i]){
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkFirstBackRerouting(i);
					} else {
						a = checkFirstBack(i);
					}
				}
				bool isSame = true;
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						if (bp_rr_prev[link[b][c]] != path_back_rr[link[b][c]][i])
						{
							isSame = false;
						}
					}
				}
				if (isSame == false)
				{
					rerouteOp++;
				}
				if(a != bp_ind[i]){
					realOp++;
					if(bpState[i]){		// Actual primary
						bpState[i] = lpState[i];
						lpState[i] = !bpState[i];
						togOp++;
					}
					mov_time++;
				}
				asignBpRerouting(i, a);
			}
		}
		ret_time += double(mov_time)*DEFRAG_TIME;
		nextEvent = eventQueue.top();
		if(nextEvent.type == 0 && (time_slot_now+ret_time >= nextEvent.time_slot_now || ret_time >= DEFRAG_TOTAL_TIME_MAX)){
			time_slot_now += ret_time;
			return 0;	
		}
		index1++;
	}
	return 0;
}


void addToList(int index, int lpsort, int bpsort) {
	lpNode *newlpNode = new lpNode;
	lpNode *newbpNode = new lpNode;
	lpNode *newlpNode2 = new lpNode;
	lpNode *newbpNode2 = new lpNode;

	newlpNode->x = index;							//lp index
	newlpNode->y = lpsort;							//lp length or size or both, for sorting purpose
	newbpNode->y = bpsort;
	newbpNode->x = index;
	newlpNode->next = NULL;
	newbpNode->next = NULL;

	lpNode **list;
	lpNode **bplist;
	lpNode **listall;

	list = &activeList;
	bplist = &backupList;
	listall = &mixtList;

	if(*list == NULL || (*list)->y < newlpNode->y ) {
		newlpNode->next = *list;
		*list = newlpNode;
	}
	else{
		lpNode *cur = *list;				// The Cur points to the first lpNode
		while(cur->next != NULL && cur->next->y <= newlpNode->y ) {
			cur = cur->next;
		}
		newlpNode->next = cur->next;
		cur->next = newlpNode;
	}

	if(*bplist == NULL || (*bplist)->y < newbpNode->y ) {
		newbpNode->next = *bplist;
		*bplist = newbpNode;
	}
	else{
		lpNode *cur = *bplist;				// The Cur points to the first lpNode
		while(cur->next != NULL && cur->next->y <= newbpNode->y ) {
			cur = cur->next;
		}
		newbpNode->next = cur->next;
		cur->next = newbpNode;
	}

	// Adding both primary and backup to the same list(listall)
	newlpNode2->x = index;
	newlpNode2->y = lpsort;
	newbpNode2->y = bpsort;
	newbpNode2->x = index;
	newlpNode2->next = NULL;
	newbpNode2->next = NULL;
	newlpNode2->z = 1;
	newbpNode2->z = 0;

	if(*listall == NULL || (*listall)->y < newlpNode2->y ) {
		newlpNode2->next = *listall;
		*listall = newlpNode2;
		newbpNode2->next = newlpNode2->next;
		newlpNode2->next = newbpNode2;
	}
	else{
		lpNode *cur = *listall;				// The Cur points to the first lpNode
		while(cur->next != NULL && cur->next->y <= newlpNode2->y ) {
			cur = cur->next;
		}
		newlpNode2->next = cur->next;
		cur->next = newlpNode2;

		lpNode *cur2 = *listall;				// The Cur points to the first lpNode
		while(cur2->next != NULL && cur2->next->y <= newbpNode2->y ) {
			cur2 = cur2->next;
		}
		newbpNode2->next = cur2->next;
		cur2->next = newbpNode2;
	}
	return;
}

void addToList2(int listId, int index, int st){

	lpNode *newlpNode = new lpNode;//メモリの確保 intみたいな感じ

	newlpNode->x = index;							//lp index
	newlpNode->z = st;
	newlpNode->next = NULL;

	lpNode **list;//*listはポインタ, **listはポインタのポインタ

	if(listId == 1) list = &realList;//アドレスを代入
	if(listId == 2) list = &tempList;

	if(*list == NULL){//もしrealListまたはtempListがなければ
		newlpNode->next = *list;//nextポインタをNULLにする
		*list = newlpNode;//realListやtempListのポインタを新しく定める
	}
	else{//もしrealListまたはtempListがあれば
		lpNode *cur = *list;// The Cur points to the first lpNode
		while(cur->next) cur = cur->next;//realListまたはtempList構造体の末尾に移動
		cur->next = newlpNode;//realListまたはtempList構造体の末尾に追加
	}
	return;
}

void delFromList(int a, int n)
{
	lpNode *currP, *prevP;
	lpNode **list;
 	
 	if(a == 1) list = &activeList;
	if(a == 2) list = &mixtList;
	if(a == 3) list = &realList;

  prevP = NULL;
  for (currP = *list; currP != NULL; prevP = currP, currP = currP->next) {
    if (currP->x == n) {
      if (prevP == NULL) {
        *list = currP->next;
      } else {
        prevP->next = currP->next;
      }

      free(currP);

      break;
    }
  }

   for (currP = *list; currP != NULL; prevP = currP, currP = currP->next) {
    if (currP->x == n) {
      if (prevP == NULL) {
        *list = currP->next;
      } else {
        prevP->next = currP->next;
      }
      free(currP);
      return;
    }
  }
}

void delFromList2(int a, int n, int st)
{
  lpNode *currP, *prevP;

  lpNode **list;
  if(a == 1) list = &activeList;
	if(a == 2) list = &mixtList;
	if(a == 3) list = &realList;
	if(a == 4) list = &tempList;

  /* For 1st node, indicate there is no previous. */
  prevP = NULL;

  for (currP = *list; currP != NULL; prevP = currP, currP = currP->next) { // Protection path
		 //リストの末尾に達するまで
	if (currP->x == n && currP->z ==st) {  /* Found it. *///リストの中にあったならば
	  if (prevP == NULL) {//Listの初っ端にあったならば
		*list = currP->next;//Listの最初の位置をずらす
	  } else {
		prevP->next = currP->next;//前の構造体にとっての次のポインタを, 現在の構造体のポインタにする.
	  }
	  free(currP);//現在の構造体のメモリを解放する
	  return;
	}
  }

  return;
}

int firstFit1_1(int lp)
{
	int a, b;
	a = checkFirstPrimRerouting(lp);
	if(a == INF){
		blocked++;
		return 0;
	}
	b = checkFirstBackRerouting(lp);
	if(b == INF){
		blocked++;
		return 0;
	}
	asignRerouting(lp, a);
	asignBpRerouting(lp, b);
	isactive[lp] = 1;
	return 1;
}

int removeLP1_1(int lp)
{
	int s = source[lp], d= dest[lp];
	int index = spec_ind[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];

	if(a){
		for(i=0; i<b; i++){
			for(j=0;j<LINK_NUM;j++){
				if (path_prim_rr[j][lp] == 1 && spec[index+i][j] == 0)
				{
					cout << "path_prim_rr[" << j << "][" << lp << "] = " << path_prim_rr[j][lp] << ", spec[" << index+i << "][" << j << "] = " << spec[index+i][j] << endl;
					throw "[error] cannot delete primary path";
				}
				spec[index+i][j] = path_prim_rr[j][lp] ^ spec[index+i][j];
			}
		}
		index = bp_ind[lp];
		if(index < INF){
			for(i=0; i<b; i++){
				for(j=0;j<LINK_NUM;j++){
					if (path_back_rr[j][lp] == 1 && spec[index+i][j] == 0)
					{
						cout << "path_back_rr[" << j << "][" << lp << "] = " << path_back_rr[j][lp] << ", spec[" << index+i << "][" << j << "] = " << spec[index+i][j] << endl;
						throw "[error] cannot delete backup path";
					}
					spec[index+i][j] = path_back_rr[j][lp] ^ spec[index+i][j];
				}
			}
		}
		isactive[lp] = 0;
	}
	return 0;
}

int checkFirstBack(int lp)
{
	bool asigned = 0, nofit = 0;
	int b= lp_size[lp];
	int index = 0;

	while(index <= (CAPASITY-b) && !asigned)
	{
		if (!b) break;
		int i;
		for (i = 0; i < LINK_NUM; i++)
		{
			if (spec[index][i] && path_back_rr[i][lp])
			{
				index++;
				break;
			} else {
				nofit = 0;
				for (int j = 0; j < b; j++)
				{
					if (spec[index+j][i] && path_back_rr[i][lp])
					{
						index += j;
						nofit = 1;
						break;
					}
				}
			}
			if (nofit) break;
		}
		if (i == LINK_NUM && !nofit){
			asigned= 1;
			return index;
		}
	}
	if(!asigned)
		return INF;
	return 0;
}

int checkFirstBackRerouting(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;
	bool isGetRoot=0;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = 0;

	while(index <= (CAPASITY-b) && !asigned)   		  //Checking all spectrum range
	{
		if(!b) break;
		isGetRoot = searchBackRoute(index, lp);
		if(isGetRoot){
			asigned= 1;
			return index;
		}else{
			index++;
		}
	}
	if(!asigned)
		return INF;
	return 0;
}

int asignBp(int lp, int index)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	bp_ind[lp] = index;
	for(j=0;j<b;j++){
		for(p=0;p<LINK_NUM;p++){
			if(spec[index+j][p] == 1 && path_back[s][d][p] ==1){
				cout << "spec[" <<  index+j << "][" << p << "] = " << spec[index+j][p] << endl;
				cout << "path_back[" << s << "][" << d << "][" << p << "] = " << path_back[s][d][p] << endl;
				throw "バックアップパス割り当てエラー";
			}
			spec[index+j][p] = spec[index+j][p] || path_back[s][d][p];
		}
	}
	return 0;
}

int asignBpRerouting(int lp, int index)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];

	bp_ind[lp] = index;
	for(j=0;j<b;j++){
		for(p=0;p<LINK_NUM;p++){
			if(spec[index+j][p] == 1 && path_back_rr[p][lp] ==1) throw "バックアップパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || path_back_rr[p][lp];
		}
	}

	return 0;
}

int initialize(void)
{
	for (int i = 0; i < NODE_NUM; i++)
	{
		for (int j = 0; j < NODE_NUM; j++)
		{
			for (int p = 0; p < LINK_NUM; p++)
			{
				path_prim[i][j][p] = 0;
				path_back[i][j][p] = 0;
			}
			hops[i][j] = 0;
			part[i][j] = 0;
			link[i][j] = INF;
			for (int k = 0; k < NODE_NUM; k++)
			{
				for (int l = 0; l < NODE_NUM; l++)
				{
					linked_path[i][j][k][l] 		= 0;
					linked_bp[i][j][k][l] 			= 0;
					linked_crosspath[i][j][k][l]	= 0;
				}
			}
		}
	}

	for (int i = 0; i < LINK_NUM; i++)
	{
		for (int j = 0; j < REQ_NUM; j++)
		{
			path_prim_rr[i][j] 	= 0;
			path_back_rr[i][j] 	= 0;
		}
	}

	return 0;
}

int reInitialize(void)
{
	blocked=0;
	togOp =0;
	realOp = 0;
	rerouteOp = 0;

	for(int i=0;i<REQ_NUM;i++){
		spec_ind[i]=0; isactive[i]=0;
		lpState[i]=1; bpState[i]=0;
	}

	// make priority empty
	while(!eventQueue.empty()){
		eventQueue.pop();
	}
	while(!deleteQueue.empty()){
		deleteQueue.pop();
	}

	for(int i=0;i<CAPASITY;i++){
		for(int j=0;j<LINK_NUM;j++)  spec[i][j]= 0;
	}

	activeList = NULL;
	backupList = NULL;
	mixtList = NULL;

	return 0;
}

int initializeEvent(void)
{
	for (int i = 0; i < REQ_NUM; i++)
	{
		endEvent[i].time_slot_now = t_exp_event[i];
		endEvent[i].type = 1;
		endEvent[i].lpNum = i;
		startEvent[i].time_slot_now = t_req_event[i];
		startEvent[i].type = 0;
		startEvent[i].lpNum = i;
	}

	defragEvent.clear();

	defragCount = round(endEvent[REQ_NUM-1].time_slot_now / DEFRAG_INTERVAL);
	defragEvent.resize(defragCount);
	for (int j = 0; j < defragCount ; j++)
	{
		defragEvent[j].time_slot_now = j * DEFRAG_INTERVAL;
		defragEvent[j].type = 2;
	}

	return 0;
}

int genDemands(int load)
{
	double expired_num = 1/double(HOLDING_TIME);
	double inter_arr = double(HOLDING_TIME)/load;

	srand (seed2);

	default_random_engine generator (seed1);
	exponential_distribution<double> next_arr(1/inter_arr);
	exponential_distribution<double> hold_time_gen(expired_num);
	uniform_int_distribution<int> traff_dist(1, REQ_SIZE_MAX);

	ofstream ofs_input;
    ofs_input.open("./../result/input.txt");
	if(!ofs_input){
		cout << "[error] cannot open input file"<< endl;
		return 1;
	}

	ofs_input << "expired number per 1 sec and inter arrival time_slot_now:=" << expired_num <<", "<<  inter_arr << endl;
	ofs_input << "lp number, source, destination, size, arrival time_slot_now, holding time_slot_now:=" << endl;

	t_req_event[0] = 0;
	for (int i = 0; i < REQ_NUM; ++i)
	{
		lp_size[i] = traff_dist(generator);
		source[i] = rand() % NODE_NUM;
		dest[i] = rand() % NODE_NUM;
		while (source[i] == dest[i]) dest[i] = rand() % NODE_NUM;
		t_hold_event[i] = hold_time_gen(generator);
		t_exp_event[i] = t_req_event[i] + t_hold_event[i];
		t_req_event[i+1] = t_req_event[i] + next_arr(generator);
		ofs_input << i  << ": "<< source[i]  <<" " << dest[i]  <<" " << lp_size[i]  <<" " << t_req_event[i] <<" " << t_hold_event[i] << endl;
		for (int j = 0; j< NODE_NUM; j++)
		{
			for (int k = 0; k < NODE_NUM; k++)
			{
				if (link[j][k] < LINK_NUM){
					path_prim_rr[link[j][k]][i] = path_prim[source[i]][dest[i]][link[j][k]];
					path_back_rr[link[j][k]][i] = path_back[source[i]][dest[i]][link[j][k]];
				}
			}
		}
	}

	ofs_input << ":;"<< endl << endl;
	ofs_input.close();
	return 0;
}

int checkFirstPrim(int lp)
{
	bool asigned = 0, nofit = 0;
	int b = lp_size[lp];
	int index = 0;

	while (index <= (CAPASITY-b) && !asigned)
	{
		if (!b) break;
		int i;
		for (i = 0; i < LINK_NUM; i++)
		{
			if (spec[index][i] && path_prim_rr[i][lp]){
				index++;
				break;
			} else {
				nofit = 0;
				for(int j = 0; j < b; j++)
				{
					if(spec[index+j][i] && path_prim_rr[i][lp]){
						index += j;
						nofit = 1;
						break;
					}
				}
			}
			if (nofit) break;
		}
		if (i == LINK_NUM && !nofit){
			asigned= 1;
			return index;
		}
	}
	if(!asigned)
		return INF;
	return 0;
}

int checkFirstPrimRerouting(int lp)
{
	bool asigned = 0;
	bool isGetRoot = 0;
	int b= lp_size[lp];
	int index=0;

	while(index <= (CAPASITY-b) && !asigned)
	{
		if (!b) break;
		isGetRoot = searchPrimRoute(index, lp);
		if(isGetRoot){
			asigned= 1;
			return index;
		} else {
			index++;
		}
	}
	if(!asigned) return INF;
	return 0;
}


int printSpec()
{
	int i,j;
	cout << "1:low index " << CAPASITY - 1 << endl;
	cout << " Spectrum :=" << endl;
	cout << "  l :";
	for (i=0;i<LINK_NUM;i++) cout <<"  "<< i ;
	cout << endl;
	for (i=CAPASITY-1; i>=0; i--){
		if(i / 10 < 1) cout << " ";
		if(i / 100 <1) cout << " ";
		cout << i << " :";
		for(j=0;j<LINK_NUM;j++){
			cout << "  " << spec[i][j];
		}
		cout << endl;
	}
	cout << "  l :";
	for (i=0;i<LINK_NUM;i++) cout <<"  "<< i ;
	cout << endl;
	cout << endl;
	return 0;
}

int asign(int lp, int index)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];

	spec_ind[lp] = index;

	for(j=0;j<b;j++){						//Asigning
		for(p=0;p<LINK_NUM;p++){
			if(spec[index+j][p] == 1 && path_prim[s][d][p] ==1){
				cout << "index + j  = " << index + j << endl;
				cout << "link_num = " << p << endl;
				throw "プライマリパス割り当てエラー";
			}
			spec[index+j][p] = spec[index+j][p] || path_prim[s][d][p];
		}
	}
	return 0;

}

int asignRerouting(int lp, int index)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];

	spec_ind[lp] = index;
	for(j=0;j<b;j++){						//Asigning
		for(p=0;p<LINK_NUM;p++){
			if(spec[index+j][p] == 1 && path_prim_rr[p][lp] ==1) throw "プライマリパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || path_prim_rr[p][lp];
		}
	}

	return 0;
}

int checkExactPrim(int lp)
{
	int b= lp_size[lp];
	int nonalign = 0, consecSB =0;

	for(int i=0;i<CAPASITY;i++){
		nonalign = 0;
		for(int j=0;j<LINK_NUM;j++){
			if(path_prim_rr[j][lp]){
				if(spec[i][j]) nonalign =1;
			}
		}
		if(nonalign==0) consecSB++;
		if(nonalign){
			if(consecSB == b) return i-b;
			consecSB = 0;
		}
		if(i==CAPASITY-1){
			if(consecSB == b) return i-b+1;
			consecSB = 0;
		}
	}
	return CAPASITY;
}

int checkExactBack(int lp)
{
	int b = lp_size[lp];
	int nonalign = 0, consecSB =0;

	for(int i = 0; i < CAPASITY; i++)
	{
		nonalign = 0;
		for(int j = 0; j < LINK_NUM; j++)
		{
			if (path_back_rr[j][lp]){
				if(spec[i][j]) nonalign =1;
			}
		}
		if(nonalign==0) consecSB++;
		if(nonalign)
		{
			if(consecSB == b) return i-b;
			consecSB = 0;
		}
		if(i==CAPASITY-1){
			if(consecSB == b) return i-b+1;
			consecSB = 0;
		}
	}
	return CAPASITY;
}

int checkExactPrimRerouting(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;
	bool isGetRoot;

	for(i=0;i<CAPASITY-b;i++){
		isGetRoot = 0;
		isGetRoot = searchPrimRoute(i, lp);
		if(isGetRoot){
			isGetRoot = 0;
			isGetRoot = searchPrimRoute(i+1, lp);
			if(!isGetRoot){
				// cout << "return i = " << i << endl;
				return i;
			}else{
				i++;//無駄に探さない
			}
		}
	}
	return CAPASITY;
}

int checkExactBackRerouting(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;
	bool isGetRoot;

	for(i=0;i<CAPASITY-b;i++){
		isGetRoot = 0;
		isGetRoot = searchBackRoute(i, lp);
		if(isGetRoot){
			isGetRoot = 0;
			isGetRoot = searchBackRoute(i+1, lp);
			if(!isGetRoot){
				// cout << "return i = " << i << endl;
				return i;
			}else{
				i++;//無駄に探さない
			}
		}
	}
	return CAPASITY;
}

int searchPrimRoute(int s, int lp)
{
	int i, j, k;
	bool unconnectedFlag;
	int a, b;
	int to, cost;

	priority_queue<Node, vector<Node>, greater<Node> > q; // 優先度付き待ち行列

	Node Nodes[NODE_NUM];
	Node doneNode;

	//ノード情報input
	for(i=0; i<NODE_NUM; i++){
		for(j=0; j<NODE_NUM; j++){
			if(link[i][j]>LINK_NUM) continue;
			unconnectedFlag = 0;
			for(k=0; k<lp_size[lp]; k++){
				if(spec[s+k][link[i][j]] == 1 || path_back_rr[link[i][j]][lp]){//使われているまたはバックアップパスが使っている
					unconnectedFlag = 1;
				}
			}
			if(!unconnectedFlag){
				Nodes[i].edges_to.push_back(j);
				Nodes[i].edges_cost.push_back(1);

			}
		}
	}

	// 初期化
	for(i=0; i<NODE_NUM; i++){
		Nodes[i].done = false;
		Nodes[i].cost = -1;
		Nodes[i].nodeNum = i;
		Nodes[i].from = INF;
	}

	Nodes[source[lp]].cost = 0; // スタートノードのコストは0

	q.push(Nodes[source[lp]]);

	// アルゴリズム実行
	while(!q.empty() && Nodes[dest[lp]].from == INF){
		// 確定ノードを取り出す
		doneNode = q.top();
		q.pop();
		//既に確定されているか確認
		if(doneNode.done) continue;
	  	// 確定フラグを立てる
	  	doneNode.done = true;
	  	// 接続先ノードの情報を更新する
	  	for(i = 0; i < doneNode.edges_to.size(); i++){
	    	to = doneNode.edges_to[i];
	    	cost = doneNode.cost + doneNode.edges_cost[i];
			    if(Nodes[to].cost < 0 || cost < Nodes[to].cost){
			    	Nodes[to].cost = cost;
					Nodes[to].from = doneNode.nodeNum;//接続元ノードを示す
					q.push(Nodes[to]);
				}
		}
	}


	//ルートの代入
	a = dest[lp];
	b = Nodes[a].from;

	// limit hop number
	int hop_counter = 0;
	while (b < NODE_NUM && a != b) {
		a = b;
		b = Nodes[a].from;
		hop_counter++;
	}
	if (hop_counter > MAX_HOP_NUM)
	{
		return 0;
	}

	//ルートの代入
	a = dest[lp];
	b = Nodes[a].from;

	if (b < NODE_NUM) {
		for (j = 0; j < NODE_NUM; j++) {//initialize
			for (k = 0; k < NODE_NUM; k++) {
				if(link[j][k]>LINK_NUM)continue;
				path_prim_rr[link[j][k]][lp] = 0;
			}
		}
		// bool isSame = true;
		while (b < NODE_NUM && a != b) {
			path_prim_rr[link[b][a]][lp] = 1;
			// std::cout << "path_prim_rr[" << link[b][a] << "][" << lp << "] = " << path_prim_rr[link[b][a]][lp] << '\n';
			a = b;
			b = Nodes[a].from;
		}
		return 1;//割り当て確定
	}else{
		return 0;//割り当て不可
	}

}

int searchBackRoute(int s, int lp)
{
	int i, j, k;
	bool unconnectedFlag;
	int a, b;
	int to, cost;

	priority_queue<Node, vector<Node>, greater<Node> > q; // 優先度付き待ち行列

	Node Nodes[NODE_NUM];
	Node doneNode;

	//ノード情報input
	for(i=0; i<NODE_NUM; i++){
		for(j=0; j<NODE_NUM; j++){
			if(link[i][j]>LINK_NUM) continue;
			unconnectedFlag = 0;
			for(k=0; k<lp_size[lp]; k++){
				if(spec[s+k][link[i][j]] == 1 || path_prim_rr[link[i][j]][lp]){//使われているまたはプライマリパスが使っている
					// cout << "path_prim_rr[" << link[i][j] << "][" << lp<< "] = " << path_prim_rr[link[i][j]][lp] << endl;
					unconnectedFlag = 1;
				}
			}
			if(!unconnectedFlag){
				Nodes[i].edges_to.push_back(j);
				Nodes[i].edges_cost.push_back(1);
				// Nodes[j].edges_to.push_back(i);
				// Nodes[j].edges_cost.push_back(1);
			}
		}
	}

	for(i= 0; i < NODE_NUM ; i++){
	    for(j=0;j < Nodes[i].edges_to.size() ; j++){
	      // cout << "Nodes[" << i << "] = " << Nodes[i].edges_to[j] << '\n';
	    }
 	}

	// 初期化
	for(i=0; i<NODE_NUM; i++){
		Nodes[i].done = false;
		Nodes[i].cost = -1;
		Nodes[i].nodeNum = i;
		Nodes[i].from = INF;
	}

	Nodes[source[lp]].cost = 0; // スタートノードのコストは0

	q.push(Nodes[source[lp]]);

	// アルゴリズム実行
	while(!q.empty() && Nodes[dest[lp]].from == INF){
	  	// 確定ノードを取り出す
		doneNode = q.top();
	 	q.pop();
		//既に確定されているか確認
		if(doneNode.done) continue;
	  	// 確定フラグを立てる
	  	doneNode.done = true;
	  	// 接続先ノードの情報を更新する
	  	for(i = 0; i < doneNode.edges_to.size(); i++){
	    	to = doneNode.edges_to[i];
	    	cost = doneNode.cost + doneNode.edges_cost[i];
	    	if(Nodes[to].cost < 0 || cost < Nodes[to].cost){
	    		Nodes[to].cost = cost;
				Nodes[to].from = doneNode.nodeNum;//接続元ノードを示す
				q.push(Nodes[to]);
			}
		}
	}

	//ルートの代入
	a = dest[lp];
	b = Nodes[a].from;

	// limit hop number
	int hop_counter = 0;
	while (b < NODE_NUM && a != b) {
		a = b;
		b = Nodes[a].from;
		hop_counter++;
	}
	if (hop_counter > MAX_HOP_NUM)
	{
		return 0;
	}

	//ルートの代入
	a = dest[lp];
	b = Nodes[a].from;

	if (b < NODE_NUM && a != b) {
		for (j = 0; j < NODE_NUM; j++) {//initialize
			for (k = 0; k < NODE_NUM; k++) {
				if(link[j][k]>LINK_NUM)continue;
				path_back_rr[link[j][k]][lp] = 0;
			}
		}
		while (b < NODE_NUM && a != b) {
			path_back_rr[link[b][a]][lp] = 1;
			// std::cout << "path_back_rr[" << link[b][a] << "][" << lp << "] = " << path_back_rr[link[b][a]][lp] << '\n';
			a = b;
			b = Nodes[a].from;
		}
		return 1;//割り当て確定
	}else{
		return 0;//割り当て不可
	}
}

void deleteLPRerouting(int lp, int p)
{
	int s = source[lp], d= dest[lp];
	int index = spec_ind[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];


	if(a){
		if(p==0 || p==1){
			for(i=0; i<b; i++){
				for(j=0;j<LINK_NUM;j++){
					if(spec[index+i][j] == 0 && path_prim_rr[j][lp] == 1) throw "プライマリパス消去エラー";
					spec[index+i][j] = path_prim_rr[j][lp] ^	spec[index+i][j];
				}
			}
		}

		if(p==0 || p==2){
			index = bp_ind[lp];
			if(index == INF) return;
			for(i=0; i<b; i++){
				for(j=0;j<LINK_NUM;j++){
					if(spec[index+i][j] == 0 && path_back_rr[j][lp] == 1) throw "バックアップパス消去エラー";
					spec[index+i][j] = path_back_rr[j][lp] ^ spec[index+i][j];
				}
			}
		}
	}
}


