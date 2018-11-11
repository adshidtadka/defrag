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
int checkFirstFit(int);
int checkFirstFitRerouting(int);
int retuneDownNonRerouting_0();
int retuneDownRerouting_0();
int asign(int, int);
int asignRerouting(int,int);
int asignBp(int, int);
int asignBpRerouting(int, int);
int printSpec(void);
int genDemands(int);
void delFromList(int, int);
void addToList(int, int, int);
int checkExactFit(int);
int checkExactBp(int);
int checkExactFitRerouting(int);
int checkExactBpRerouting(int);
int getPrimRoot(int, int);
int getBackRoot(int, int);
void deleteLP(int, int);
void deleteLPRerouting(int, int);
int removeLP1_1(int, int);
int firstFit1_1(int, int);
int checkFirstBp(int);
int checkFirstBpRerouting(int);
int retuneBp(int);
int readResultPy();
int readResultReroutingPy();
void statDefragPy(int);
void statDefragReroutingPy(int);
int statAlgo(void);
int statReroutingAlgo(void);
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
bool spec[CAPASITY][LINK_NUM], path[NODE_NUM][NODE_NUM][LINK_NUM];
bool linked_path[NODE_NUM][NODE_NUM][NODE_NUM][NODE_NUM];
bool linked_bp[NODE_NUM][NODE_NUM][NODE_NUM][NODE_NUM];
bool linked_crosspath[NODE_NUM][NODE_NUM][NODE_NUM][NODE_NUM];
int blocked;
int isactive[REQ_NUM];
int t_req[REQ_NUM], t_hold[REQ_NUM], t_exp[REQ_NUM];
double t_req_event[REQ_NUM], t_hold_event[REQ_NUM], t_exp_event[REQ_NUM];
int hops[NODE_NUM][NODE_NUM], bhops[NODE_NUM][NODE_NUM];
int link[NODE_NUM][NODE_NUM];
int last_lp;
int midlp = REQ_NUM;
unsigned seed1 = 123;
unsigned seed2 =  125;
struct lpNode{
	int x;
	int y;
	int z;
	lpNode *next;
};
lpNode *activeList = NULL;
lpNode *midFitList = NULL;
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
bool path_rr[LINK_NUM][REQ_NUM], bp_rr[LINK_NUM][REQ_NUM];
int part[NODE_NUM][NODE_NUM];
int algoCall;
int t_temp =0;
int last_ret= 0;
int temp_max = DEFRAG_TOTAL_TIME_MAX;
double t;
int bp_ind[REQ_NUM];
bool bp[NODE_NUM][NODE_NUM][LINK_NUM];
int togOp;
int realOp;
int rerouteOp;
bool lpState[REQ_NUM];
bool bpState[REQ_NUM];
int last_blocked = 0;
struct Event
{
	double time;
	int lpNum;
	int type;	
	bool operator> (const Event &event) const {
   	return (time > event.time);
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
	int i,j,k,l, p;
	int a, b;
	char tmp= 'r';

	ifstream fin;

	//	Register network
	stringstream ss;
	string s;
	ss << "./../network/net_" << argv[1] << ".txt";
	s = ss.str();
	fin.open (s);//ノードの接続関係が EM:=で与えられている
	if (!fin){
		cout <<"Cannot open network model" << endl;
		return 1;
	}

	fin.ignore(INT_MAX,'=');
	for(k=0; k<LINK_NUM; k++){
		fin >> a >> b ;
		link[a][b]=k;
	}
	fin.close();

//	Register routing table primary
	ss.str("");
	ss << "./../network/dp" << argv[1] << "_x_result1.txt";
	s = ss.str();
	fin.open (s); //行数をカウントしてlに代入する
	if (!fin)
	{
		cout <<"Cannot open routing file" << endl;
		return 1;
	}

	fin.ignore(INT_MAX,'=');//=が出るまで無視
	//	fin.ignore(INT_MAX,':');
	fin >> tmp;					//tmp1 is a char

	l=0;
	while(tmp !=';'){//もし最終行でなければ
		l++;
		fin.ignore(INT_MAX,'\n');//改行まで無視
		fin >> tmp; //次の行の最初の文字を読み取る
	}
	fin.close();

	ss.str("");
	ss << "./../network/dp" << argv[1] << "_x_result1.txt";
	s = ss.str();
	fin.open (s); //今度こそプライマリのパス情報を得る path[i][j][p] hops[i][j]
	if (!fin)
	{
		cout <<"Cannot open routing file" << endl;
		return 1;
	}

	//cout << l <<": lines" << endl;
	fin.ignore(INT_MAX,'='); //=まで無視
	for (k=0;k<l;k++){ //行数分forループ
		fin >> i >> j >> a >> b;
		fin.ignore(INT_MAX,'\n'); //改行まで無視する
		p= link[a][b]; //pにリンク番号を代入
		// cout << p <<": "<< a << b << endl;
		path[i][j][p] = 1; /*path[i][j][p]は発ノードと着ノードを結ぶパスが使う
		リンク番号を代入すると1をとるバイナリ変数*/
	//	if(k<20) cout << i <<" "<< j <<" " << p << endl;　
		++hops[i][j];  			// Counting number of hops
		//hops[i][j]発ノードと着ノードを結ぶパスが何ホップでできているか
	}
	fin.close();

//	Register routing table backup
	ss.str("");
	ss << "./../network/dp" << argv[1] << "_x_result2.txt";
	s = ss.str();
	fin.open (s);//行数をカウントしてlに代入してる
	{
		if (!fin)
		{
			cout <<"Cannot open routing file" << endl;
			return 1;
		}

		fin.ignore(INT_MAX,'=');
	//	fin.ignore(INT_MAX,':');
		fin >> tmp;					//tmp1 is a char

		l=0;
		while(tmp !=';'){
			l++;
			fin.ignore(INT_MAX,'\n');
			fin >> tmp;
		}
	}
	fin.close();

	ss.str("");
	ss << "./../network/dp" << argv[1] << "_x_result2.txt";
	s = ss.str();
	fin.open (s);//今度こそプライマリのパス情報を得る bp[i][j][p] bhops[i][j]
	{
		if (!fin)
		{
			cout <<"Cannot open routing file" << endl;
			return 1;
		}

		//cout << l <<": lines" << endl;
		fin.ignore(INT_MAX,'=');
		for (k=0;k<l;k++){
			fin >> i >> j >> a >> b;
			fin.ignore(INT_MAX,'\n');
			p= link[a][b];
		//	cout << p <<": "<< a << b << endl;
			bp[i][j][p] = 1;
			// if(k<20) cout << i <<" "<< j <<" " << p << endl;
			++bhops[i][j];  			// Counting number of hops
			//bhops[i][j]発ノードと着ノードを結ぶバックアップパスパスが何ホップでできているか
		}
	}
	fin.close();

	for (i=0;i<NODE_NUM;i++){				//Comparing path, may be usefull
		for (j=0;j<NODE_NUM;j++){
			for (k=0;k<NODE_NUM;k++){
				for (l=0;l<NODE_NUM;l++){
					if(k == i && l == j) linked_path[i][j][k][l] = 0; //2組の発着ノードの組が同じ発着ノードであれば無視する
					else{
						for(p=0;p<LINK_NUM;p++){ //リンクの数だけforループを回す
							if(path[i][j][p] && path[k][l][p]) linked_path[i][j][k][l] = 1;
							//linked_path[i][j][k][l]はプライマリパスに関して2組の発着ノードを代入したとき
							//同じリンクを使っているようであれば1をとるバイナリ変数
							if(bp[i][j][p] && bp[k][l][p]) linked_bp[i][j][k][l] = 1; //
							//linked_bp[i][j][k][l]はバックアップパスに関して2組の発着ノードを代入したとき
							//同じリンクを使っているようであれば1をとるバイナリ変数
							if(path[i][j][p] && bp[k][l][p]) linked_crosspath[i][j][k][l] = 1;
							//linked_crosspath[i][j][k][l]はプライマリパスの発着ノードとバックアップパスの発着ノードを代入したとき
							//同じリンクを使っているようであれば1をとるバイナリ変数
						}
					}
				}
			}
		}
	}
	genDemands(load);
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
			if(path[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
		for(j=0;j<LINK_NUM;j++){
			if(bp[s][d][j]) ofs2 << ind <<" "<< j << endl;
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
			if(path[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
		for(j=0;j<LINK_NUM;j++){
			if(bp[s][d][j]) ofs2 << ind <<" "<< j << endl;
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
					if (path_rr[link[i][j]][lp]){
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
					if (bp_rr[link[i][j]][lp]){
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
	if(!algoCall) return 0;

	if(algoCall==1){
		statAlgo();
		retuneDownNonRerouting_0();
		return 0;
	}

	if(algoCall==2){
		statReroutingAlgo();
		retuneDownRerouting_0();
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
			// printSpec();
			asign(lp, b);
	//	if(a<6 || a> 54) cout << "LP " << lp << " at " << b << endl;
			fin.ignore(INT_MAX,'\n');
			fin >> a >> b ;
			asignBp(lp, b);
	//	if(a<6 || a> 54) cout << "LP " << lp << " at " << b << endl;
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
				path_rr[i][lp] = 0;
				bp_rr[i][lp]   = 0;
			}

			fin.ignore(INT_MAX,'=');
			fin >> a >> b >> c >> d;
			while(!d){//同じパスなら
				path_rr[link[b][c]][lp] = 1;
				fin.ignore(INT_MAX,'\n');
				// cout << "path_rr[link[" << b << "][" << c << "]][" << lp << "] = " << path_rr[link[b][c]][lp] << endl;
				fin >> a >> b >> c >> d;
			}
			fin.ignore(INT_MAX,'=');
			fin >> a >> b >> c >> d;
			while(!d){//同じパスなら
				bp_rr[link[b][c]][lp] = 1;
				fin.ignore(INT_MAX,'\n');
				// cout << "bp_rr[link[" << b << "][" << c << "]][" << lp << "] = " << bp_rr[link[b][c]][lp] << endl;
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

int statAlgo()
{
	int a, b, lp, lp2;
	double ret_time = 0;
	int realMov = 1;
	int s1, d1, s2, d2;
	int st;

	lpNode *cur;
	lpNode *cur1;
	lpNode *cur2;

	tempList = NULL;
	cur = mixtList;
	while ( cur != NULL ) {//mixtListの内容をtempListの末尾に追加
		addToList2(2, cur->x, cur->z);//1ならrealList, 2ならtempList
	//	cout << "cur x, z= " << cur->x << ", " << cur->z << endl;
		cur= cur->next;
	}
	cur = tempList;

	//デフラグ終了時刻を出す
	double fin_time;
	fin_time = finTime();
	if (fin_time < 0) // queue is empty
	{
		return 0;
	}

	while(realMov){//最初はrealMov=1 パスの割り当てがあれば
		realMov = 0;
		tempList = NULL;
		cur = mixtList;
		while ( cur != NULL ) { //mixtListの内容をtempListの末尾に追加
			addToList2(2, cur->x, cur->z);//1ならrealList, 2ならtempList
	//		cout << "cur x, z= " << cur->x << ", " << cur->z << endl;
			cur= cur->next;
		}
		cur = tempList;
		realList = NULL;
		cur1 = NULL;
	//	cout << "cur = tempList " << endl<< endl;
		while ( cur != NULL ) {	// Checking all active LPs　tempListがなくなるまで
			realList = NULL;
			cur1 = cur;//cur1は現在のtempList
			cur = cur->next;//curが次のtempListを指す
		//	cout << "cur1 = cur " << endl<< endl;

			while ( cur1 != NULL ) { //tempListがなくなるまで
				//cout << "cur1 = cur 2" << endl<< endl;
				lp = cur1->x;//現在のtempListの中身
				st = cur1->z;//現在のtempListの中身
				cur1 = cur1->next;//cur1が次のtempListを指す

				cur2 = realList;
				//cout << "cur2 = realList , cur1->x, cur1->z " << lp <<", " << st << endl;
				//if(cur2) cout << "cur2 = realList 16, cur2->x " << realList->x << endl;

				if(cur2 == NULL){//realListがなければ
					//	cout << "cur2 = realList 101 , reallist->x " << realList->x << endl<< endl;
					addToList2(1, lp, st);//realListの末尾に現在のtempListの中身を追加
					if(cur && cur->x == lp && cur->z == st) cur = cur->next;
					//tempListに次の構造体が存在していて内容が1つ前と一致していたならば, 次の次の構造体へ
					delFromList2(4, lp, st);//追加したlpをtempListから削除
					//		cout << "cur2 = realList 102, reallist->x " << realList->x << endl;
				}else{//realListがあれば
					//	cout << "cur2 = realList 102, reallist->x " << realList->x << endl;
					//	cout << "cur2 = realList 11 " << endl<< endl;
					int conf = 0;
					while(cur2 != NULL){//realListがtempListのリンクとかぶっていないか
					//	cout << "cur2 = realList 12, reallist->x " << realList->x << endl;
						s2 = source[cur2->x];//realListの発ノード
						d2= dest[cur2->x];//realListの着ノード
						cur2 = cur2->next;//次のrealListへ
						//cout << "cur2 = realList 15, cur2->x " << realList->x << endl;
						s1 = source[lp];//tempListの発ノード
						d1= dest[lp];//tempListの着ノード
						//cout << "cur2 = realList 13, cur2->x " << endl;
						if(linked_path[s1][d1][s2][d2] || linked_bp[s1][d1][s2][d2] || linked_crosspath[s1][d1][s2][d2]){
							//もし2組の発着ノードに関して同じリンクを使っているならば
							conf =1;
							//cout << " conf " << endl;
							break;//かぶっているとしてチェック終了
						}
					}
					//cout << "cur2 = realList 14 " << realList->x << endl<< endl;
					if(conf==0){//もしrealListがtempListのリンクとかぶっていなれば
						addToList2(1, lp, st);//realListの末尾に現在のtempListの中身を追加
						if(cur && cur->x == lp && cur->z == st) cur = cur->next;
						//tempListに次の構造体が存在していて内容が1つ前と一致していたならば, 次の次の構造体へ
						delFromList2(4, lp, st);//追加したlpをtempListから削除
						//1:active 2:mixtList 3:realList 4:tempList
					}
					//cout << "cur2 = realList 14 " << endl<< endl;
				}//realListがあれば
			}//tempListがなくなるまで

			//cout << "cur2 = realList 2" << endl<< endl;
			cur2 = realList;
			int realcheck = realOp;//帯域移動操作の総数が変わっているかあとで確認
			while ( cur2 != NULL ){//realListがなくなるまで
				lp = cur2->x;
				b = cur2->z;
				cur2 = cur2->next;
				//cout << "cur2 = realList 1, lp= " << lp << endl<< endl;
				if(b){//プライマリパスならば
					//cout << "cur2 = realList 2, lp= " << lp << endl<< endl;
					deleteLP(lp, 1);// Remove LP from spec to avoid self-conflict
					//プライマリパスのみ消去
					// printSpec();
					//a = checkFirstFit(lp);				// Return spec_ind[lp] at worst
					a = checkExactFit(lp);
					if(a==CAPASITY || a > spec_ind[lp]) a = checkFirstFit(lp);
					//もしスロット番号が小さくならないようであれば
					//if(a > spec_ind[lp]) a = spec_ind[lp];
					if(a != spec_ind[lp]){//スロット番号が小さくなっていれば(INFならブロッキング)
						realOp++;//帯域移動操作の総数
						realMov++;
						if(lpState[lp]){		// Actual primary
							bpState[lp] = lpState[lp];
							lpState[lp] = !bpState[lp];
							togOp++;
						}
					}
					asign(lp, a);
				}
				//cout << "cur2 = realList 3" << endl<< endl;
				if(!b){//バックアップパスならば
					deleteLP(lp, 2);            // Remove LP from spec to avoid self-conflict
					//バックアップパスのみ消去
					//printSpec();
					// a = checkFirstBp(lp);				// Return spec_ind[lp] at worst
					a = checkExactBp(lp);
					if(a==CAPASITY || a > spec_ind[lp]) a = checkFirstBp(lp);
					if(a != bp_ind[lp]){//スロット番号が小さくなっていれば(INFならブロッキング)
						realOp++;
						realMov++;
						if(bpState[lp]){		// Actual primary
							bpState[lp] = lpState[lp];
							lpState[lp] = !bpState[lp];
							togOp++;
						}
					}
					asignBp(lp, a);
				}
			}//realListがなくなる
			cur2 = realList;//realListの先頭に戻す
			//cout << "cur2 = realList 3" << endl<< endl;
			while ( cur2 != NULL ){//realListがなくなるまで
				lp = cur2->x;
				cur2 = cur2->next;
				delFromList(3, lp);//1:active 2:mixtList 3:realList 4:tempList
			}			
			
			if(realcheck != realOp){//帯域移動操作数が変わっていれば
				ret_time += DEFRAG_TIME;					//increment defrag time
				if((t+ret_time >= fin_time || ret_time >= temp_max) || eventQueue.empty()){
					// t += ret_time;cout << "t = " << t << ", ret_time = " << ret_time << ", fin_time = " << fin_time << endl;	
					return 0;	
				}
				//新しいパスがきているか, 100000ステップを超えたならばデフラグ終了
			}
		}//tempListがなくなるまで
	}//realMovが0以外なら続くのwhile文 帯域移動がなければデフラグ終了
	// t += ret_time;
	return 0;
}

int statReroutingAlgo()
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
	while ( cur != NULL ) {//mixtListの内容をtempListの末尾に追加
		addToList2(2, cur->x, cur->z);//1ならrealList, 2ならtempList
	//	cout << "cur x, z= " << cur->x << ", " << cur->z << endl;
		cur= cur->next;
	}
	cur = tempList;

	//デフラグ終了時刻を出す
	double fin_time;
	fin_time = finTime();
	if (fin_time < 0) // queue is empty
	{
		return 0;
	}
	// cout << "t = " << t << ", ret_time = " << ret_time << ", fin_time = " << fin_time << endl;

	while(realMov){//最初はrealMov=1 パスの割り当てがあれば
		realMov = 0;
		tempList = NULL;
		cur = mixtList;
		while ( cur != NULL ) { //mixtListの内容をtempListの末尾に追加
			addToList2(2, cur->x, cur->z);//1ならrealList, 2ならtempList
	//		cout << "cur x, z= " << cur->x << ", " << cur->z << endl;
			cur= cur->next;
		}
		cur = tempList;
		realList = NULL;
		cur1 = NULL;
	//	cout << "cur = tempList " << endl<< endl;
		while ( cur != NULL ){	// Checking all active LPs　tempListがなくなるまで
			realList = NULL;
			cur1 = cur;//cur1は現在のtempList
			cur = cur->next;//curが次のtempListを指す
			//	cout << "cur1 = cur " << endl<< endl;

			while ( cur1 != NULL ) { //tempListがなくなるまで
			//		cout << "cur1 = cur 2" << endl<< endl;
				lp = cur1->x;//現在のtempListの中身
				st = cur1->z;//現在のtempListの中身
				cur1 = cur1->next;//cur1が次のtempListを指す

				cur2 = realList;
			//		cout << "cur2 = realList , cur1->x, cur1->z " << lp <<", " << st << endl;
			//	if(cur2) cout << "cur2 = realList 16, cur2->x " << realList->x << endl;

				if(cur2 == NULL){//realListがなければ
					// cout << "cur2 = realList 101 , reallist->x " << realList->x << endl<< endl;
					addToList2(1, lp, st);//realListの末尾に現在のtempListの中身を追加
					if(cur && cur->x == lp && cur->z == st) cur = cur->next;
					//tempListに次の構造体が存在していて内容が1つ前と一致していたならば, 次の次の構造体へ
					delFromList2(4, lp, st);//追加したlpをtempListから削除
					// cout << "cur2 = realList 102, reallist->x " << realList->x << endl;
				}else{//realListがあれば
					//	cout << "cur2 = realList 102, reallist->x " << realList->x << endl;
					//	cout << "cur2 = realList 11 " << endl<< endl;
					int conf = 0;
					while(cur2 != NULL){//realListがなくなるまで、tempListのリンクとかぶっていないか確認
					//	cout << "cur2 = realList 12, reallist->x " << realList->x << endl;

						if(lpState[cur2->x]){//プライマリ
							if(lpState[lp]){
								for(i=0;i<LINK_NUM;i++){
									if(path_rr[i][cur2->x] && path_rr[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<LINK_NUM;i++){
									if(path_rr[i][cur2->x] && bp_rr[i][lp])conf =1;
								}

							}
						}else{//バックアップ
							if(lpState[lp]){
								for(i=0;i<LINK_NUM;i++){
									if(bp_rr[i][cur2->x] && path_rr[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<LINK_NUM;i++){
									if(bp_rr[i][cur2->x] && bp_rr[i][lp]) conf =1;
								}
							}
						}
						if(conf) break;//かぶっているとしてチェック終了
						cur2 = cur2->next;//次のrealListへ
					}//realListがなくなるまで
				//	cout << "cur2 = realList 14 " << realList->x << endl<< endl;
					if(conf==0){//もしrealListがtempListのリンクとかぶっていなれば
						addToList2(1, lp, st);//realListの末尾に現在のtempListの中身を追加
						if(cur && cur->x == lp && cur->z == st) cur = cur->next;
						//tempListに次の構造体が存在していて内容が1つ前と一致していたならば, 次の次の構造体へ
						delFromList2(4, lp, st);//追加したlpをtempListから削除
						//1:active 2:mixtList 3:realList 4:tempList
					}
					//	cout << "cur2 = realList 14 " << endl<< endl;
				}//realListがあれば
			}//tempListがなくなるまで

			// cout << "cur2 = realList 2" << endl<< endl;
			cur2 = realList;
			int realcheck = realOp;//帯域移動操作の総数が変わっているかあとで確認
			while ( cur2 != NULL ){//realListがなくなるまで
				lp = cur2->x;
				b =  lpState[lp];
				cur2 = cur2->next;
				if(b){//プライマリパスならば
					int path_rr_prev[LINK_NUM]; //to compare the previous root and new route
					for (int i = 0; i < LINK_NUM; ++i)
					{
						path_rr_prev[i] = 0;
					}
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							path_rr_prev[link[i][j]] = path_rr[link[i][j]][lp];
						}
					}
					deleteLPRerouting(lp, 1);// Remove LP from spec to avoid self-conflict
					//プライマリパスのみ消去
					a = checkExactFitRerouting(lp);
					// std::cout << " after checkExactFitRerouting a = " << a << '\n';
					if(a==CAPASITY || a > spec_ind[lp]){
						a = checkFirstFitRerouting(lp);
						// std::cout << " after checkFirstFitRerouting a = " << a << '\n';
					}
					bool isSame = true; //check the rerouted or not rerouted
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							if (path_rr_prev[link[i][j]] != path_rr[link[i][j]][lp])
							{
								isSame = false;
							}
						}
					}
					if (isSame == false)
					{
						rerouteOp++;
					}
					if(a != spec_ind[lp]){//スロット番号が小さくなっていれば(INFならブロッキング)
						realOp++;//帯域移動操作の総数
						realMov++;
						if(lpState[lp]){		// Actual primary
							bpState[lp] = lpState[lp];
							lpState[lp] = !bpState[lp];
							togOp++;
						}
					}
					asignRerouting(lp, a);
				}
			//	cout << "cur2 = realList 3" << endl<< endl;
				if(!b){//バックアップパスならば
					int bp_rr_prev[LINK_NUM]; //to compare the previous root and new route
					for (int i = 0; i < LINK_NUM; ++i)
					{
						bp_rr_prev[i] = 0;
					}
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							bp_rr_prev[link[i][j]] = bp_rr[link[i][j]][lp];
						}
					}
					// std::cout << "バックアップ, lp = " << lp << ", bp_ind[lp] = " << bp_ind[lp] << ", source[lp] = " << source[lp] << ", dest[lp] =" << dest[lp] << ", lp_size[lp] =" <<lp_size[lp]<<'\n';
					deleteLPRerouting(lp, 2);            // Remove LP from spec to avoid self-conflict
					//バックアップパスのみ消去
					a = checkExactBpRerouting(lp);
					// std::cout << " after checkExactFitRerouting a = " << a << '\n';
					if(a==CAPASITY || a > bp_ind[lp]){
						a = checkFirstBpRerouting(lp);
						// std::cout << " after checkFirstFitRerouting a = " << a << '\n';
					}
					bool isSame = true; //check the rerouted or not rerouted
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							if (bp_rr_prev[link[i][j]] != bp_rr[link[i][j]][lp])
							{
								isSame = false;
							}
						}
					}
					if (isSame == false)
					{
						rerouteOp++;
					}
					// std::cout << "after checkFirstBpRerouting a = " << a << '\n';
					if(a != bp_ind[lp]){//スロット番号が小さくなっていれば(INFならブロッキング)
						realOp++;
						realMov++;
						if(bpState[lp]){		// Actual primary
							bpState[lp] = lpState[lp];
							lpState[lp] = !bpState[lp];
							togOp++;
						}
					}
					asignBpRerouting(lp, a);
				}
			}//realListがなくなる
			cur2 = realList;//realListの先頭に戻す
			//	cout << "cur2 = realList 3" << endl<< endl;
			while ( cur2 != NULL ){//realListがなくなるまで
				lp = cur2->x;
				cur2 = cur2->next;
				delFromList(3, lp);//1:active 2:mixtList 3:realList 4:tempList
			}
			if(realcheck != realOp){//帯域移動操作数が変わっていれば
				ret_time += DEFRAG_TIME;					//increment defrag time
				if((t+ret_time >= fin_time || ret_time >= temp_max) || eventQueue.empty()){
					// t += ret_time;
					return 0;	
				}
				//新しいパスがきているか, 100000ステップを超えたならばデフラグ終了
			}

		}//tempListがなくなるまで
	}//realMovが0以外なら続くのwhile文 帯域移動がなければデフラグ終了
	return 0;
}

double finTime(){

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

	return eventQueue.top().time;
}

int retuneDownNonRerouting_0()
{
	int s, d;
	int index1;
	int i,j,p;
	int k=0;
	int a;

	double ret_time = 0 ;
	int mov_time = 0 ;
//	int ret_time = t ;

	index1 = 1;

	while(index1 < CAPASITY-1 && k++ < INF){			//Sweep the spectrum
		lpNode *cur = activeList; 
		mov_time = 0 ;
		while ( cur != NULL ) {						// Checking all active LPs
			i = cur->x;
			cur = cur->next;
			if(isactive[i] == 1 && spec_ind[i] == index1){
				deleteLP(i, 1);// Remove LP from spec to avoid self-conflict
				//プライマリパスのみ消去
				//a = checkFirstFit(lp);				// Return spec_ind[lp] at worst
				a = checkExactFit(i);
				if(a==CAPASITY || a > spec_ind[i]) a = checkFirstFit(i);
				//もしスロット番号が小さくならないようであれば
				//if(a > spec_ind[lp]) a = spec_ind[lp];
				if(a != spec_ind[i]){//スロット番号が小さくなっていれば(INFならブロッキング)
					realOp++;//帯域移動操作の総数
					if(lpState[i]){		// Actual primary
						bpState[i] = lpState[i];
						lpState[i] = !bpState[i];
						togOp++;
					}
					mov_time++;
				}
				asign(i, a);
			}

			if(isactive[i] == 1 && bp_ind[i] == index1){
				deleteLP(i, 2);            // Remove LP from spec to avoid self-conflict
				//バックアップパスのみ消去
				a = checkFirstBp(i);				// Return spec_ind[lp] at worst
				if(a != bp_ind[i]){//スロット番号が小さくなっていれば(INFならブロッキング)
					realOp++;
					if(bpState[i]){		// Actual primary
						bpState[i] = lpState[i];
						lpState[i] = !bpState[i];
						togOp++;
					}
					mov_time++;
				}
					asignBp(i, a);
			}
		}
		//increment defrag time
		ret_time += double(mov_time)*DEFRAG_TIME;
		nextEvent = eventQueue.top(); 	// next event
		if(nextEvent.type == 0 && (t+ret_time >= nextEvent.time || ret_time >= temp_max)){
			// t += ret_time;
			return 0;	
		}
		index1++;
	}
	return 0;
}

int retuneDownRerouting_0()
{
	int s, d;
	int index1;
	int i,j,p;
	int k=0;
	int a;

	int ret_time = 0 ;
	int mov_time = 0 ;
//	int ret_time = t ;

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
						path_rr_prev[link[b][c]] = path_rr[link[b][c]][i];
					}
				}
				deleteLPRerouting(i, 1);// Remove LP from spec to avoid self-conflict
				//プライマリパスのみ消去
				a = checkExactFitRerouting(i);
				// std::cout << " after checkExactFitRerouting a = " << a << '\n';
				if(a==CAPASITY || a > spec_ind[i]){
					a = checkFirstFitRerouting(i);
					// std::cout << " after checkFirstFitRerouting a = " << a << '\n';
				}
				bool isSame = true; //check the rerouted or not rerouted
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						if (path_rr_prev[link[b][c]] != path_rr[link[b][c]][i])
						{
							isSame = false;
						}
					}
				}
				if (isSame == false)
				{
					rerouteOp++;
				}
				if(a != spec_ind[i]){//スロット番号が小さくなっていれば(INFならブロッキング)
					realOp++;//帯域移動操作の総数
					if(lpState[i]){		// Actual primary
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
						bp_rr_prev[link[b][c]] = bp_rr[link[b][c]][i];
					}
				}
				deleteLPRerouting(i, 2);            // Remove LP from spec to avoid self-conflict
				//バックアップパスのみ消去
				a = checkExactBpRerouting(i);
				// std::cout << " after checkExactFitRerouting a = " << a << '\n';
				if(a==CAPASITY || a > bp_ind[i]){
					a = checkFirstBpRerouting(i);
				}
				bool isSame = true; //check the rerouted or not rerouted
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						if (bp_rr_prev[link[b][c]] != bp_rr[link[b][c]][i])
						{
							isSame = false;
						}
					}
				}
				if (isSame == false)
				{
					rerouteOp++;
				}
				if(a != bp_ind[i]){//スロット番号が小さくなっていれば(INFならブロッキング)
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
		//increment defrag time
		ret_time += double(mov_time)*DEFRAG_TIME;
		nextEvent = eventQueue.top(); 	// next event
		if(nextEvent.type == 0 && (t+ret_time >= nextEvent.time || ret_time >= temp_max)){
			t += ret_time;
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

//	cout <<"2: At t= " << t << " Allocating demand " << index <<endl;

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
//	cout <<"3: At t= " << t << " Allocating demand " << index <<endl;

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
//	cout <<"4: At t= " << t << " Allocating demand " << index <<endl;
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

//	cout <<"2: At t= " << t << " Allocating demand " << index <<endl;

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
  lpNode *currP, *prevP; //2つのポインタ宣言

  lpNode **list;//ポインタのポインタ
  if(a == 1) list = &activeList; //activeListのアドレスをlistに代入
	if(a == 2) list = &mixtList;	//
	if(a == 3) list = &realList;
	//どのリストを変更するかをaでもらう

  /* For 1st node, indicate there is no previous. */
  prevP = NULL;

  /*
   * Visit each node, maintaining a pointer to
   * the previous node we just visited.
   */
  for (currP = *list; currP != NULL; prevP = currP, currP = currP->next) {
	//listがなくなるまで
    if (currP->x == n) {  /* Found it. lp番号が一致 */
      if (prevP == NULL) {
        /* Fix beginning pointer. */
        *list = currP->next;
      } else {
        /*
         * Fix previous node's next to
         * skip over the removed node.
         */
        prevP->next = currP->next;
      }

      /* Deallocate the node. */
      free(currP);

      /* Done searching. */
      break;
    }
  }

   for (currP = *list; currP != NULL; prevP = currP, currP = currP->next) {
    if (currP->x == n) {  /* Found it. */
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

int firstFit1_1(int lp, int algoCall)
{
	int a, b;
	if(algoCall == 2 || algoCall == 4){
		a = checkFirstFitRerouting(lp);					// For primary

		if(a == INF){
			blocked++;
			return 0;
		}

		b = checkFirstBpRerouting(lp);					// For Backup
		if(b == INF){
			blocked++;
			return 0;
		}
		//ブロッキングが起こらなかったらパスを割り当てる
		asignRerouting(lp, a);
		asignBpRerouting(lp, b);
		isactive[lp] = 1;

		return 1;
	}else{
		a = checkFirstFit(lp);					// For primary

		if(a == INF){
			blocked++;
			return 0;
		}

		b = checkFirstBp(lp);					// For Backup
		if(b == INF){
			blocked++;
			return 0;
		}
		//ブロッキングが起こらなかったらパスを割り当てる
		asign(lp, a);
		asignBp(lp, b);
		isactive[lp] = 1;

		return 1;
	}
}

int removeLP1_1(int lp, int algoCall) //切断するパスのlpindexをもらう
{
	int s = source[lp], d= dest[lp];
	int index = spec_ind[lp],  b= lp_size[lp];//占有する帯域スロット数
	int i,j;
	int a = isactive[lp];

//	cout << "loop0, a=" << a << endl;
	if(algoCall == 2 || algoCall == 4){
		if(a){										// a xor a =0, therefore if a lp is active
			for(i=0; i<b; i++){//占有する帯域スロットの回数繰り返す									// an xor with its path will remove it
				for(j=0;j<LINK_NUM;j++){//全てのリンクについて
					spec[index+i][j] = path_rr[j][lp] ^ spec[index+i][j];
					//全てのリンクの占有する帯域スロット番号について
					//プライマリパスが通っているところはビット反転させる
				}
			}

			index = bp_ind[lp];
			if(index < INF){
				for(i=0; i<b; i++){									// an xor with its path will remove it
					for(j=0;j<LINK_NUM;j++){
						spec[index+i][j] = bp_rr[j][lp] ^ spec[index+i][j];
						//全てのリンクの占有する帯域スロット番号について
						//バックアップパスが通っているところはビット反転させる
					}
				}
			}
			isactive[lp] = 0;//非アクティブにする
		}
		return 0;

	}else{
		if(a){										// a xor a =0, therefore if a lp is active
			for(i=0; i<b; i++){//占有する帯域スロットの回数繰り返す									// an xor with its path will remove it
				for(j=0;j<LINK_NUM;j++){//全てのリンクについて
					spec[index+i][j] = path[s][d][j] ^ spec[index+i][j];
					//全てのリンクの占有する帯域スロット番号について
					//プライマリパスが通っているところはビット反転させる
				}
			}

			index = bp_ind[lp];
			if(index < INF){
				for(i=0; i<b; i++){									// an xor with its path will remove it
					for(j=0;j<LINK_NUM;j++){
						spec[index+i][j] = bp[s][d][j] ^ spec[index+i][j];
						//全てのリンクの占有する帯域スロット番号について
						//バックアップパスが通っているところはビット反転させる
					}
				}
			}
			isactive[lp] = 0;//非アクティブにする
		}
		return 0;
	}
}

int checkFirstBp(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = 0;

	while(index <= (CAPASITY-b) && !asigned)   		  //Checking all spectrum range
	{
		if(!b) break;

		for(i=0;i<LINK_NUM;i++){							//Check path availibility
			if(spec[index][i] && bp[s][d][i]){
				index++;
				break;
			}else{
				nofit = 0;
				for(j=0;j<b;j++){							//Checking if it fit (size)
					if(spec[index+j][i] && bp[s][d][i]){
						index += j;
						nofit = 1;
						break;
					}
				}
			}
			if(nofit) break;
		}
		if(i==LINK_NUM && !nofit){
			asigned= 1;
			return index;
		}
	}
	if(!asigned)
		return INF;
	return 0;
}

int checkFirstBpRerouting(int lp)
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
		isGetRoot = getBackRoot(index, lp);
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
			if(spec[index+j][p] == 1 && bp[s][d][p] ==1){
				cout << "spec[" <<  index+j << "][" << p << "] = " << spec[index+j][p] << endl;
				cout << "bp[" << s << "][" << d << "][" << p << "] = " << bp[s][d][p] << endl;
				throw "バックアップパス割り当てエラー";
			}
			spec[index+j][p] = spec[index+j][p] || bp[s][d][p];
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
			if(spec[index+j][p] == 1 && bp_rr[p][lp] ==1) throw "バックアップパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || bp_rr[p][lp];
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
				path[i][j][p] = 0;
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
			path_rr[i][j] 	= 0;
			bp_rr[i][j] 	= 0;
		}
	}

	return 0;
}

int reInitialize(void)
{
	int i,j;
	blocked=0;
	togOp =0;
	realOp = 0;
	rerouteOp = 0;

	for(i=0;i<REQ_NUM;i++){
		spec_ind[i]=0; isactive[i]=0;
		lpState[i]=1; bpState[i]=0;
	}

	//make priority empty
	while(!eventQueue.empty()){
		eventQueue.pop();
	}
	while(!deleteQueue.empty()){
		deleteQueue.pop();
	}

	for(i=0;i<CAPASITY;i++){
		for(j=0;j<LINK_NUM;j++)  spec[i][j]= 0;
	}

	activeList = NULL;
	midFitList = NULL;
	backupList = NULL;
	mixtList = NULL;

	return 0;
}

int initializeEvent(void)						// Set everything to zero save routings and demands
{
	int i,j;

	for(i=0;i<REQ_NUM;i++){
		endEvent[i].time  = 0;
		endEvent[i].type  = 0;
		endEvent[i].lpNum = 0;
		startEvent[i].time  = 0;
		startEvent[i].type  = 0;
		startEvent[i].lpNum = 0;
	}

	defragEvent.clear();

	return 0;
}

int genDemands(int load)
{
	int i,j,k;
	int arr_int;
	double arr_int_event;
	double temp1;

	double mu = 1/double(HOLDING_TIME);
	double inter_arr = double(HOLDING_TIME)/load;

	srand (seed2);

	default_random_engine generator (seed1);
	exponential_distribution<double> next_arr(1/inter_arr);
	exponential_distribution<double> hold_time(mu);
	uniform_int_distribution<int> traff_dist(1, REQ_SIZE_MAX);

	ofstream ofs_result_txt;
    ofs_result_txt.open("./../result/input_demands1.txt");
	if(!ofs_result_txt){
		cout<< "Cannot open input_demands1 file"<<endl;
		return 1;
	}

	ofstream ofs2;
    ofs2.open("./../result/input_demands2.txt");
	if(!ofs2){
		cout<< "Cannot open input_demands2 file"<<endl;
		return 1;
	}

	ofs_result_txt << "mu and inter:=" << mu <<", "<<  inter_arr<< endl;
	ofs_result_txt << "LP number, s, d, size, Arrival time, holding time:=" << endl;
	ofs2 << "mu and inter:=" << mu <<", "<<  inter_arr << endl;
	ofs2 << "LP number, s, d, size, Arrival time, holding time:=" << endl;

	int max_hold = 0;
	t_req[0]= 0;
	t_req_event[0] = 0;

	for (int i=0; i<REQ_NUM; ++i){ //リクエストの数だけ
		temp1 = hold_time(generator); //指数関数分布で継続時間を乱数として得る
		t_hold[i] = int(temp1)+1;		//100をかけて1を足したものが継続時間となる
		t_hold_event[i] = temp1;
		if(max_hold < t_hold[i]){
			max_hold = t_hold[i];//継続時間の最大値を超えていれば最大値とする
		}
		arr_int_event = next_arr(generator);
		arr_int = int(arr_int_event)+1; //ポアソン分布で到着間隔を得る
		lp_size[i] = traff_dist(generator); //等確率で占有帯域スロット数を得る
		source[i] = rand() %NODE_NUM;//発ノードをランダムで得る
		dest[i] = rand() %NODE_NUM;//着ノードをランダムに得る
			while(source[i] == dest[i]) dest[i] = rand() %NODE_NUM;//発ノードと着ノードが同じにならないようにする
		t_exp[i] = t_req[i] + t_hold[i]; //継続時間と到着時刻を合わせて切断時刻を計算
		t_exp_event[i] = t_req_event[i] + t_hold_event[i]; //event driven
		ofs_result_txt << i  << ": "<< source[i]  <<" " << dest[i]  <<" " << lp_size[i]  <<" " << t_req[i] <<" " << t_hold[i] << endl;
		ofs2 << i  << ": "<< source[i]  <<" " << dest[i]  <<" " << lp_size[i]  <<" " << t_req_event[i] <<" " << t_hold_event[i] << endl;
		t_req[i+1]= t_req[i] + arr_int; //到着時刻と到着間隔を合わせて次の到着時刻を計算
		t_req_event[i+1] = t_req_event[i] + arr_int_event; //event driven
		//lp番号ごとのルートを持つ変数も定める
		for(j=0;j<NODE_NUM;j++){
			for(k=0;k<NODE_NUM;k++){
				if(link[j][k] < LINK_NUM){
					path_rr[link[j][k]][i] = path[source[i]][dest[i]][link[j][k]];
					bp_rr[link[j][k]][i] = bp[source[i]][dest[i]][link[j][k]];
				}
			}
		}
	}

	ofs_result_txt << ":;"<< endl << endl;
	ofs_result_txt.close();
	ofs2 << ":;"<< endl << endl;
	ofs2.close();
	return 0;
}

int checkFirstFit(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = 0;

	while(index <= (CAPASITY-b) && !asigned) //Checking all spectrum range 全ての帯域スロットを調べる
	{
		if(!b) break;//bが0ならブレイク

	//	cout << "LP checking index: " << lp << ", "<< index << endl;

		for(i=0;i<LINK_NUM;i++){		//Check path availibility
			if(spec[index][i] && path[s][d][i]){ //リンクが使われていたならば
				index++;
				break;
			}else{//リンクが空いていたならば
				nofit = 0;
				for(j=0;j<b;j++){	//Checking if it fit (size)
					if(spec[index+j][i] && path[s][d][i]){ //リンクが使われていたならば
						index += j;//無駄にさがさない
						nofit = 1;
						break;
					}
				}
			}
			if(nofit) break;
		}
		if(i==LINK_NUM && !nofit){//もし全てのリンクを調べ終え、リンクが使われていなければ
			asigned= 1;//割り当て確定
		//	cout << "LP index " << lp << " will be first-asigned to index " << index << endl;
		//	cout << "LP index is allocated #" << isactive[lp] << endl;
			return index;
		}
	}
	if(!asigned)//もし割り当てされていなければ
		return INF;//ブロッキングが起こったということ
	return 0;
}

int checkFirstFitRerouting(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;
	bool isGetRoot = 0;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = 0;

	while(index <= (CAPASITY-b) && !asigned) //Checking all spectrum range 全ての帯域スロットを調べる
	{
		if(!b) break;//bが0ならブレイク
		//	cout << "LP checking index: " << lp << ", "<< index << endl;
		isGetRoot = getPrimRoot(index, lp);
		//	cout << "LP checked index: " << lp << ", "<< index << endl;
		if(isGetRoot){//もし全てのリンクを調べ終え、リンクが使われていなければ
			asigned= 1;//割り当て確定
		//	cout << "LP index " << lp << " will be first-asigned to index " << index << endl;
		//	cout << "LP index is allocated #" << isactive[lp] << endl;
			return index;
		}else{
			index++;
		}
	}
	if(!asigned)//もし割り当てされていなければ
		return INF;//ブロッキングが起こったということ
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
			if(spec[index+j][p] == 1 && path[s][d][p] ==1){
				cout << "index + j  = " << index + j << endl;
				cout << "link_num = " << p << endl;
				throw "プライマリパス割り当てエラー";
			}
			spec[index+j][p] = spec[index+j][p] || path[s][d][p];
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
			if(spec[index+j][p] == 1 && path_rr[p][lp] ==1) throw "プライマリパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || path_rr[p][lp];
		}
	}

	return 0;
}

int checkExactFit(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;

	double specFrag =0;
	double totalfrag=0, avefrag=0;
	int totalhops=0;

	int maxSB = 0, avSB = 0;
	int nonalign = 0, consecSB =0, maxcons=0;

	for(i=0;i<CAPASITY;i++){		// Checking spectrum for available aligned SB
		nonalign =0 ;
		for(j=0;j<LINK_NUM;j++){//全てのリンクについて
			if(path[s][d][j]){		// Path s-d using link j
				if(spec[i][j])  nonalign =1;	// SB i not aligned through the path of s-d
			}
		}//もしあるリンクが使用されていたならばnonalign =1となる
		if(nonalign==0) consecSB++;//Increasing consecutive aligned SB
		if(nonalign){			// End of consecutive aligned SB exactfitしていたら
			if(consecSB == b)//もし必要な帯域スロット数と同じならば
				return i-b;
			consecSB = 0;
		}
		if(i==CAPASITY-1){//もし最高帯域スロットに達したならば
			if(consecSB == b)//もし必要な帯域スロットと同じならば
				return i-b+1;
			consecSB = 0;
		}
	}

	return CAPASITY;
}

int checkExactBp(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;

	double specFrag =0;
	double totalfrag=0, avefrag=0;
	int totalhops=0;

	int maxSB = 0, avSB = 0;
	int nonalign = 0, consecSB =0, maxcons=0;

	for(i=0;i<CAPASITY;i++){		// Checking spectrum for available aligned SB
		nonalign =0 ;
		for(j=0;j<LINK_NUM;j++){//全てのリンクについて
			if(bp[s][d][j]){		// Path s-d using link j
				if(spec[i][j])  nonalign =1;	// SB i not aligned through the path of s-d
			}
		}//もしあるリンクが使用されていたならばnonalign =1となる
		if(nonalign==0) consecSB++;//Increasing consecutive aligned SB
		if(nonalign){			// End of consecutive aligned SB exactfitしていたら
			if(consecSB == b)//もし必要な帯域スロット数と同じならば
				return i-b;
			consecSB = 0;
		}
		if(i==CAPASITY-1){//もし最高帯域スロットに達したならば
			if(consecSB == b)//もし必要な帯域スロットと同じならば
				return i-b+1;
			consecSB = 0;
		}
	}

	return CAPASITY;
}

int checkExactFitRerouting(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;
	bool isGetRoot;

	for(i=0;i<CAPASITY-b;i++){
		isGetRoot = 0;
		isGetRoot = getPrimRoot(i, lp);
		if(isGetRoot){
			isGetRoot = 0;
			isGetRoot = getPrimRoot(i+1, lp);
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

int checkExactBpRerouting(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;
	bool isGetRoot;

	for(i=0;i<CAPASITY-b;i++){
		isGetRoot = 0;
		isGetRoot = getBackRoot(i, lp);
		if(isGetRoot){
			isGetRoot = 0;
			isGetRoot = getBackRoot(i+1, lp);
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

int getPrimRoot(int s, int lp)
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
				if(spec[s+k][link[i][j]] == 1 || bp_rr[link[i][j]][lp]){//使われているまたはバックアップパスが使っている
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
				path_rr[link[j][k]][lp] = 0;
			}
		}
		// bool isSame = true;
		while (b < NODE_NUM && a != b) {
			path_rr[link[b][a]][lp] = 1;
			// std::cout << "path_rr[" << link[b][a] << "][" << lp << "] = " << path_rr[link[b][a]][lp] << '\n';
			a = b;
			b = Nodes[a].from;
		}
		return 1;//割り当て確定
	}else{
		return 0;//割り当て不可
	}

}

int getBackRoot(int s, int lp)
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
				if(spec[s+k][link[i][j]] == 1 || path_rr[link[i][j]][lp]){//使われているまたはプライマリパスが使っている
					// cout << "path_rr[" << link[i][j] << "][" << lp<< "] = " << path_rr[link[i][j]][lp] << endl;
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
				bp_rr[link[j][k]][lp] = 0;
			}
		}
		while (b < NODE_NUM && a != b) {
			bp_rr[link[b][a]][lp] = 1;
			// std::cout << "bp_rr[" << link[b][a] << "][" << lp << "] = " << bp_rr[link[b][a]][lp] << '\n';
			a = b;
			b = Nodes[a].from;
		}
		return 1;//割り当て確定
	}else{
		return 0;//割り当て不可
	}

}



void deleteLP(int lp, int p)	// p=0 delete both, 1 del prim and 2 del backup
{
	int s = source[lp], d= dest[lp];
	int index = spec_ind[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];

//	cout << "loop1, index=" << index << endl;

	if(a){						// a xor a =0, therefore if a lp is active an xor with its path will remove it
		if(p==0 || p==1){
		//	if(lp==181) cout << "index= " << index << endl;
			for(i=0; i<b; i++){//占有帯域スロット
				for(j=0;j<LINK_NUM;j++){//全てのリンクに関して
					if(spec[index+i][j] == 0 && path[s][d][j] == 1) throw "プライマリパス消去エラー";
					spec[index+i][j] = path[s][d][j] ^	spec[index+i][j];//pathが1ならspecは1
				}
			}
		}

		if(p==0 || p==2){
			index = bp_ind[lp];
			if(index == INF) return;//ブロッキングしている
		//	if(lp==181) cout << "Big index= " << index << endl;
			for(i=0; i<b; i++){									// an xor with its path will remove it
				for(j=0;j<LINK_NUM;j++){
					if(spec[index+i][j] == 0 && bp[s][d][j] == 1) throw "バックアップパス消去エラー";
					spec[index+i][j] = bp[s][d][j] ^ spec[index+i][j];//pathが1ならspecは1
				}
			}
		}
	}
}

void deleteLPRerouting(int lp, int p)	// p=0 delete both, 1 del prim and 2 del backup
{
	int s = source[lp], d= dest[lp];
	int index = spec_ind[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];

//	cout << "loop1, index=" << index << endl;

	if(a){						// a xor a =0, therefore if a lp is active an xor with its path will remove it
		if(p==0 || p==1){
		//	if(lp==181) cout << "index= " << index << endl;
			for(i=0; i<b; i++){//占有帯域スロット
				for(j=0;j<LINK_NUM;j++){//全てのリンクに関して
					if(spec[index+i][j] == 0 && path_rr[j][lp] == 1) throw "プライマリパス消去エラー";
					spec[index+i][j] = path_rr[j][lp] ^	spec[index+i][j];//pathが1ならspecは1
				}
			}
		}

		if(p==0 || p==2){
			index = bp_ind[lp];
			if(index == INF) return;//ブロッキングしている
		//	if(lp==181) cout << "Big index= " << index << endl;
			for(i=0; i<b; i++){									// an xor with its path will remove it
				for(j=0;j<LINK_NUM;j++){
					if(spec[index+i][j] == 0 && bp_rr[j][lp] == 1) throw "バックアップパス消去エラー";
					spec[index+i][j] = bp_rr[j][lp] ^ spec[index+i][j];//pathが1ならspecは1
				}
			}
		}
	}

}


