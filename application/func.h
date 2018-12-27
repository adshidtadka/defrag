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
#include <map>

using namespace std;

int initialize(void);
int reInitialize(void);
int initializeEvent(void);
int readInput(int, char**, int);
int retuneDown(void);
int asignPrim(int,int);
int asignBack(int, int);
int printSpec(void);
int genDemands(int);
void delFromList(int, int);
void addToList(int, int, int);
int checkExactPrimConv(int);
int checkExactBackConv(int);
int checkExactPrimProp(int, bool);
int checkExactBackProp(int, bool);
int checkFirstPrimConv(int);
int checkFirstPrimProp(int, bool);
int checkFirstBackConv(int);
int checkFirstBackProp(int, bool);
int searchRoutePrim(int, int, bool);
int searchRouteBack(int, int, bool);
void delLp(int, int);
int removeLP1_1(int);
int setupPath(int);
int startDefrag(int);
int readResultConv();
int readResultProp();
int startAlgo(void);
int writeGivenParamConv(int);
int writeGivenParamProp(int);
void addToList2(int, int, int);
void delFromList2(int, int, int);
double sum(double, int);
double ave(double, int);
double var(double, int);
double standard(double, int);
double finTime();
int isAvailablePrim(int, int, int);
int isAvailableBack(int, int, int);

bool path_prim[LINK_NUM][REQ_NUM];
bool path_back[LINK_NUM][REQ_NUM];
int link[NODE_NUM][NODE_NUM];
bool spec[CAPASITY][LINK_NUM];

int lp_size[REQ_NUM], source[REQ_NUM], dest[REQ_NUM];
int isactive[REQ_NUM];
double t_req_event[REQ_NUM], t_hold_event[REQ_NUM], t_exp_event[REQ_NUM];
int ind_prim[REQ_NUM], ind_back[REQ_NUM];
bool state_prim[REQ_NUM];
bool state_back[REQ_NUM];
int limit_hop_prim[REQ_NUM];
int limit_hop_back[REQ_NUM];
unsigned seed1 = SEED_1;
unsigned seed2 = SEED_2;
int blocked;
int algoCall;
int defragCount;
double time_slot_now;
int togOp;
int realOp;
int rerouteOp;
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
vector<Event> defragEvent;

// map<char, int> links;
// map<char, int> nodes;
// links["5node"] 	= 5;
// nodes["5node"] 	= 12;
// links["abi"] 	= 11;
// nodes["abi"] 	= 28;
// links["euro"] 	= 11;
// nodes["euro"] 	= 52;
// links["nsf"] 	= 14;
// nodes["nsf"] 	= 44;
// links["ind"] 	= 14;
// nodes["ind"] 	= 46;
// links["jap"] 	= 25;
// nodes["jap"] 	= 84;

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
	}
	fin.close();

	return 0;
}

int writeGivenParamConv(int load)
{
	int i,j,ind=0;
	int lp,n,f0,k;

	int m = 0;
	lpNode *cur = activeList;
	while ( cur != NULL ) {
		m++;
		lp = cur->x;
		cur = cur->next;
	}

	ofstream ofs2;
	ofs2.open("./../result/ssr_lno_output.txt", ios::out);
	if(!ofs2){
		cout<< "[error] cannot open ./../result/ssr_lno_output.txt"<<endl;
		return 1;
	}

	ofs2 << "param S   := " << m   << endl;
	ofs2 << "param absP:= " << m*2 << endl;
	ofs2 << "param absF:= " << CAPASITY   << endl;
	ofs2 << "param absE:= " << LINK_NUM   << endl;
	ofs2 << "param T   := " << MAX_STEP   << endl;
	ofs2 << endl;

	ofs2 << "param : s_p k_p_init n_p f_p_init := " << endl;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		n = lp_size[lp];
		k = state_prim[lp];
		if (k)
		{
			f0 = ind_prim[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			f0 = ind_back[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			ind++;
		}else{
			f0 = ind_prim[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			f0 = ind_back[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			ind++;
		}
	}
	ofs2 << ";" << endl << endl;

	ofs2 << "set P_e := " << endl;
	ind =0;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;

		for(j=0;j<LINK_NUM;j++){
			if(path_prim[j][lp]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
		for(j=0;j<LINK_NUM;j++){
			if(path_back[j][lp]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
	}
	ofs2 << ";" << endl << endl;

	ofs2.close();
	return 0;
}

int writeGivenParamProp(int load)
{
	int i,j,ind=0;
	int lp,n,f0,k;

	int m = 0;
	lpNode *cur = activeList;
	while ( cur != NULL ) {
		m++;
		lp = cur->x;
		cur = cur->next;
	}

	ofstream ofs2;
	ofs2.open("./../result/ssrr_lno_output.txt", ios::out);
	if(!ofs2){
		cout<< "Cannot open ssrr_lno_output.txt file"<<endl;
		return 1;
	}

	ofs2 << "param S   := " << m   << endl;
	ofs2 << "param absP:= " << m*2 << endl;
	ofs2 << "param absF:= " << CAPASITY   << endl;
	ofs2 << "param absE:= " << LINK_NUM   << endl;
	ofs2 << "param absV:= " << NODE_NUM   << endl;
	ofs2 << "param T   := " << MAX_STEP   << endl;
	ofs2 << endl;

	ofs2 << "param : s_p k_p_init n_p f_p_init := " << endl;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		n = lp_size[lp];
		k = state_prim[lp];
		if (k)
		{
			f0 = ind_prim[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			f0 = ind_back[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			ind++;
		}else{
			f0 = ind_prim[lp];
			ofs2 << ind << " " <<"0 "<< n <<" "<< f0 << endl;
			f0 = ind_back[lp];
			ofs2 << ind << " " <<"1 "<< n <<" "<< f0 << endl;
			ind++;
		}
	}
	ofs2 << ";" << endl << endl;

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
	ind =0;
	cur = activeList;
	while ( cur != NULL ){
		lp = cur->x;
		cur = cur->next;
		for ( i = 0; i < NODE_NUM; i++) {
			for ( j = 0; j < NODE_NUM; j++ ){
				if(link[i][j] < LINK_NUM){
					if (path_prim[link[i][j]][lp]){
						ofs2 << ind << " " << i << " " << j << endl;
					}
				}
			}
		}

		ind++;
		for ( i = 0; i < NODE_NUM; i++) {
			for ( j = 0; j < NODE_NUM; j++ ){
				if(link[i][j] < LINK_NUM){
					if (path_back[link[i][j]][lp]){
						ofs2 << ind << " " << i << " " << j << endl;
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

		ofs2 << source[lp] << " " << dest[lp] << endl;
	}
	ofs2 << ";" << endl << endl;
	ofs2.close();
	return 0;
}


int startDefrag(int load)
{
	if(algoCall == 0)
	{
		 return 0;
	}

	if(algoCall == 1 || algoCall == 2)
	{
		startAlgo();
		retuneDown();
		return 0;
	}

	if(algoCall == 3)
	{
		writeGivenParamConv(load);
		system("python ssr_lno.py");
		readResultConv();
		return 0;
	}

	if(algoCall == 4)
	{
		writeGivenParamProp(load);
		system("python ssrr_lno.py");
		readResultProp();
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
    ofs_input.open("./../result/input.txt", ios::out);
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
	}

	ofs_input << ":;"<< endl << endl;
	ofs_input.close();
	return 0;
}

int readResultConv()
{
	int i,j;
	int a,b, lp;
	ifstream fin, fin1;

	for(i=0;i<CAPASITY;i++){
		for(j=0;j<LINK_NUM;j++)  spec[i][j]= 0;
	}

	fin.open ("./../result/ssr_lno_result.txt");
		if (!fin){
			cout <<"[error] cannot read ssr_lno_result.txt file" << endl;
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
		while ( cur != NULL ) {
			lp = cur->x;
			cur = cur->next;

			fin >> a >> b ;
			asignPrim(lp, b);
			fin.ignore(INT_MAX,'\n');
			fin >> a >> b ;
			asignBack(lp, b);
			fin.ignore(INT_MAX,'\n');
		}
	fin.close();
	return 0;
}

int readResultProp()
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
		while ( cur != NULL ) {
			lp = cur->x;
			cur = cur->next;

			for(i=0;i<LINK_NUM;i++){
				path_prim[i][lp] = 0;
				path_back[i][lp]   = 0;
			}

			fin.ignore(INT_MAX,'=');
			fin >> a >> b >> c >> d;
			while(!d){
				path_prim[link[b][c]][lp] = 1;
				fin.ignore(INT_MAX,'\n');
				fin >> a >> b >> c >> d;
			}
			fin.ignore(INT_MAX,'=');
			fin >> a >> b >> c >> d;
			while(!d){
				path_back[link[b][c]][lp] = 1;
				fin.ignore(INT_MAX,'\n');
				fin >> a >> b >> c >> d;
			}
		}

		fin.ignore(INT_MAX,'=');

		cur = activeList;
		while ( cur != NULL )
		{
			lp = cur->x;
			cur = cur->next;

			fin >> a >> b ;
			asignPrim(lp, b);
			fin.ignore(INT_MAX,'\n');
			fin >> a >> b ;
			asignBack(lp, b);
			fin.ignore(INT_MAX,'\n');
		}
	fin.close();
	return 0;
}

int startAlgo()
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
						if(state_prim[cur2->x]){
							if(state_prim[lp]){
								for(i=0;i<LINK_NUM;i++){
									if(path_prim[i][cur2->x] && path_prim[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<LINK_NUM;i++){
									if(path_prim[i][cur2->x] && path_back[i][lp])conf =1;
								}

							}
						}else{
							if(state_prim[lp]){
								for(i=0;i<LINK_NUM;i++){
									if(path_back[i][cur2->x] && path_prim[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<LINK_NUM;i++){
									if(path_back[i][cur2->x] && path_back[i][lp]) conf =1;
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
				b =  state_prim[lp];
				cur2 = cur2->next;
				if(b){
					int path_prev[LINK_NUM];
					for (int i = 0; i < LINK_NUM; ++i)
					{
						path_prev[i] = 0;
					}
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							path_prev[link[i][j]] = path_prim[link[i][j]][lp];
						}
					}
					delLp(lp, 1);
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkExactPrimProp(lp, false);
					} else {
						a = checkExactPrimConv(lp);
					}
					if (a == CAPASITY || a > ind_prim[lp])
					{
						if (algoCall == 2 || algoCall == 4)
						{
							a = checkFirstPrimProp(lp, false);
						} else {
							a = checkFirstPrimConv(lp);
						}
					}
					if (a == INF)
					{
						a  = ind_prim[lp];
					}
					bool isSame = true;
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							if (path_prev[link[i][j]] != path_prim[link[i][j]][lp])
							{
								isSame = false;
							}
						}
					}
					if (isSame == false)
					{
						rerouteOp++;
					}
					if(a != ind_prim[lp]){
						realOp++;
						realMov++;
						if(state_prim[lp]){
							state_back[lp] = state_prim[lp];
							state_prim[lp] = !state_back[lp];
							togOp++;
						}
					}
					asignPrim(lp, a);
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
							bp_rr_prev[link[i][j]] = path_back[link[i][j]][lp];
						}
					}
					delLp(lp, 2);
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkExactBackProp(lp, false);
					} else {
						a = checkExactBackConv(lp);
					}
					if (a == CAPASITY || a > ind_back[lp])
					{
						if (algoCall == 2 || algoCall == 4)
						{
							a = checkFirstBackProp(lp, false);
						} else {
							a = checkFirstBackConv(lp);
						}
					}
					if (a == INF)
					{
						a = ind_back[lp];
					}
					bool isSame = true;
					for (int i = 0; i < NODE_NUM; ++i)
					{
						for (int j = 0; j < NODE_NUM; ++j)
						{
							if (link[i][j] > LINK_NUM) continue;
							if (bp_rr_prev[link[i][j]] != path_back[link[i][j]][lp])
							{
								isSame = false;
							}
						}
					}
					if (isSame == false)
					{
						rerouteOp++;
					}
					if(a != ind_back[lp]){
						realOp++;
						realMov++;
						if(state_back[lp]){
							state_back[lp] = state_prim[lp];
							state_prim[lp] = !state_back[lp];
							togOp++;
						}
					}
					asignBack(lp, a);
				}
			}
			cur2 = realList;
			while ( cur2 != NULL ){
				lp = cur2->x;
				cur2 = cur2->next;
				delFromList(3, lp);
			}
			if(realcheck != realOp){
				ret_time += PROCESSING_TIME;
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
	return eventQueue.top().time;
}

int retuneDown()
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
			if(isactive[i] == 1 && ind_prim[i] == index1){
				int path_prev[LINK_NUM]; //to compare the previous root and new route
				for (int b = 0; b < LINK_NUM; ++b)
				{
					path_prev[b] = 0;
				}
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						path_prev[link[b][c]] = path_prim[link[b][c]][i];
					}
				}
				delLp(i, 1);
				if (algoCall == 2 || algoCall == 4)
				{
					a = checkExactPrimProp(i, false);
				} else {
					a = checkExactPrimConv(i);
				}
				if(a==CAPASITY || a > ind_prim[i]){
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkFirstPrimProp(i, false);
					} else {
						a = checkFirstPrimConv(i);
					}
				}
				if (a == INF)
				{
					a = ind_prim[i];
				}
				bool isSame = true;
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						if (path_prev[link[b][c]] != path_prim[link[b][c]][i])
						{
							isSame = false;
						}
					}
				}
				if (isSame == false)
				{
					rerouteOp++;
				}
				if(a != ind_prim[i]){
					realOp++;
					if(state_prim[i]){
						state_back[i] = state_prim[i];
						state_prim[i] = !state_back[i];
						togOp++;
					}
					mov_time++;
				}
				asignPrim(i, a);
			}

			if(isactive[i] == 1 && ind_back[i] == index1){
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
						bp_rr_prev[link[b][c]] = path_back[link[b][c]][i];
					}
				}
				delLp(i, 2);
				if (algoCall == 2 || algoCall == 4)
				{
					a = checkExactBackProp(i, false);
				} else {
					a = checkExactBackConv(i);
				}
				if(a==CAPASITY || a > ind_back[i]){
					if (algoCall == 2 || algoCall == 4)
					{
						a = checkFirstBackProp(i, false);
					} else {
						a = checkFirstBackConv(i);
					}
				}
				if (a == INF)
				{
					a = ind_back[i];
				}
				bool isSame = true;
				for (int b = 0; b < NODE_NUM; ++b)
				{
					for (int c = 0; c < NODE_NUM; ++c)
					{
						if (link[b][c] > LINK_NUM) continue;
						if (bp_rr_prev[link[b][c]] != path_back[link[b][c]][i])
						{
							isSame = false;
						}
					}
				}
				if (isSame == false)
				{
					rerouteOp++;
				}
				if(a != ind_back[i]){
					realOp++;
					if(state_back[i]){		// Actual primary
						state_back[i] = state_prim[i];
						state_prim[i] = !state_back[i];
						togOp++;
					}
					mov_time++;
				}
				asignBack(i, a);
			}
		}
		ret_time += double(mov_time)*PROCESSING_TIME;
		nextEvent = eventQueue.top();
		if(nextEvent.type == 0 && (time_slot_now+ret_time >= nextEvent.time || ret_time >= DEFRAG_TOTAL_TIME_MAX)){
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

int setupPath(int lp)
{
	int a = checkFirstPrimProp(lp, true);
	if(a == INF){
		blocked++;
		return 0;
	}
	int b = checkFirstBackProp(lp, true);
	if(b == INF){
		blocked++;
		return 0;
	}
	asignPrim(lp, a);
	asignBack(lp, b);
	isactive[lp] = 1;
	return 1;
}

int removeLP1_1(int lp)
{
	int index = ind_prim[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];

	if(a){
		for(i=0; i<b; i++){
			for(j=0;j<LINK_NUM;j++){
				if (path_prim[j][lp] == 1 && spec[index+i][j] == 0)
				{
					cout << "path_prim[" << j << "][" << lp << "] = " << path_prim[j][lp] << ", spec[" << index+i << "][" << j << "] = " << spec[index+i][j] << endl;
					throw "[error] cannot delete primary path";
				}
				spec[index+i][j] = path_prim[j][lp] ^ spec[index+i][j];
			}
		}
		index = ind_back[lp];
		if(index < INF){
			for(i=0; i<b; i++){
				for(j=0;j<LINK_NUM;j++){
					if (path_back[j][lp] == 1 && spec[index+i][j] == 0)
					{
						cout << "path_back[" << j << "][" << lp << "] = " << path_back[j][lp] << ", spec[" << index+i << "][" << j << "] = " << spec[index+i][j] << endl;
						throw "[error] cannot delete backup path";
					}
					spec[index+i][j] = path_back[j][lp] ^ spec[index+i][j];
				}
			}
		}
		isactive[lp] = 0;
	}
	return 0;
}

int asignBack(int lp, int index)
{
	int i,j,p;
	int b= lp_size[lp];

	ind_back[lp] = index;
	for(j=0;j<b;j++){
		for(p=0;p<LINK_NUM;p++){
			if(spec[index+j][p] == 1 && path_back[p][lp] ==1) throw "バックアップパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || path_back[p][lp];
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
			link[i][j] = INF;
		}
	}

	for (int i = 0; i < LINK_NUM; i++)
	{
		for (int j = 0; j < REQ_NUM; j++)
		{
			path_prim[i][j] 	= 0;
			path_back[i][j] 	= 0;
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
		ind_prim[i]=0; isactive[i]=0;
		state_prim[i]=1; state_back[i]=0;
		limit_hop_prim[i] = 0; limit_hop_back[i] = 0;
	}

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
		endEvent[i].time = t_exp_event[i];
		endEvent[i].type = 1;
		endEvent[i].lpNum = i;
		startEvent[i].time = t_req_event[i];
		startEvent[i].type = 0;
		startEvent[i].lpNum = i;
	}

	defragEvent.clear();

	defragCount = round(endEvent[REQ_NUM-1].time / DEFRAG_INTERVAL);
	defragEvent.resize(defragCount);
	for (int j = 0; j < defragCount ; j++)
	{
		defragEvent[j].time = j * DEFRAG_INTERVAL;
		defragEvent[j].type = 2;
	}

	return 0;
}

int checkFirstPrimConv(int lp)
{
	for (int i = 0; i < CAPASITY - lp_size[lp]; ++i)
	{
		if (isAvailablePrim(i, 1, lp))
		{
			if (isAvailablePrim(i + 1, lp_size[lp] - 1, lp))
			{
				return i;
			}
		}
	}
	return INF;
}

int checkFirstBackConv(int lp)
{
	for (int i = 0; i < CAPASITY - lp_size[lp]; ++i)
	{
		if (isAvailableBack(i, 1, lp))
		{
			if (isAvailableBack(i + 1, lp_size[lp] - 1, lp))
			{
				return i;
			}
		}
	}
	return INF;
}

int checkFirstPrimProp(int lp, bool isSetUp)
{
	for (int i = 0; i < CAPASITY - lp_size[lp]; ++i)
	{
		if (searchRoutePrim(i, lp, isSetUp))
		{
			return i;
		}
	}
	return INF;
}

int checkFirstBackProp(int lp, bool isSetUp)
{
	for (int i = 0; i < CAPASITY - lp_size[lp]; ++i)
	{
		if (searchRouteBack(i, lp, isSetUp))
		{
			return i;
		}
	}
	return INF;
}

int asignPrim(int lp, int index)
{
	int i,j,p;
	int b= lp_size[lp];

	ind_prim[lp] = index;
	for(j=0;j<b;j++){
		for(p=0;p<LINK_NUM;p++){
			if(spec[index+j][p] == 1 && path_prim[p][lp] ==1) throw "プライマリパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || path_prim[p][lp];
		}
	}

	return 0;
}

int checkExactPrimConv(int lp)
{
	int b= lp_size[lp];
	int consecSB =0;

	for(int  i = 0; i < CAPASITY; i++)
	{
		if (isAvailablePrim(i, 1, lp))
		{
			 consecSB++;
		} 
		else 
		{
			if (consecSB == b){
				 return i-b;	
			}
			else
			{
				consecSB = 0;
			}
		}

		if (i == CAPASITY - 1)
		{
			if (consecSB == b)
			{
				 return i - b + 1;	
			}
			consecSB = 0;
		}
	}
	return CAPASITY;
}

int checkExactBackConv(int lp)
{
	int b= lp_size[lp];
	int consecSB =0;

	for(int  i = 0; i < CAPASITY; i++)
	{
		
		if (isAvailableBack(i, 1, lp))
		{
			 consecSB++;
		} 
		else 
		{
			if (consecSB == b){
				 return i-b;	
			}
			else
			{
				consecSB = 0;
			}
		}

		if (i == CAPASITY - 1)
		{
			if (consecSB == b)
			{
				 return i - b + 1;	
			}
			consecSB = 0;
		}
	}
	return CAPASITY;
}

int checkExactPrimProp(int lp, bool isSetUp)
{
	int i,j,p;
	int b= lp_size[lp];
	int lp1, index1, s1, d1;
	bool isGetRoot;

	for(i=0;i<CAPASITY-b;i++){
		isGetRoot = 0;
		isGetRoot = searchRoutePrim(i, lp, isSetUp);
		if(isGetRoot){
			// confirm non-availability
			if (i + 1 + b == CAPASITY)
			{
				return i;
			} 
			else 
			{
				if (isAvailablePrim(i + 1 + b, 1, lp) == 0)
				{
					return i;
				}
			}
		}
	}
	return CAPASITY;
}

int checkExactBackProp(int lp, bool isSetUp)
{
	int b = lp_size[lp];
	bool isGetRoot;

	for(int i = 0; i < CAPASITY - b; i++){
		isGetRoot = 0;
		isGetRoot = searchRouteBack(i, lp, isSetUp);
		if(isGetRoot){
			// confirm non-availability
			if (i + 1 + b == CAPASITY)
			{
				return i;
			} 
			else 
			{
				if (isAvailableBack(i + 1 + b, 1, lp) == 0)
				{
					return i;
				}
			}
		}
	}
	return CAPASITY;
}

int isAvailablePrim(int index, int times, int lp) {

	for (int i = index; i < index + times; ++i)
	{
		for (int j = 0; j < LINK_NUM; ++j)
		{
			if (path_prim[j][lp] && spec[i][j])
			{
				return 0;
			}
		}
	}
	return 1;
}

int isAvailableBack(int index, int times, int lp) {

	for (int i = index; i < index + times; ++i)
	{
		for (int j = 0; j < LINK_NUM; ++j)
		{
			if (path_back[j][lp] && spec[i][j])
			{
				return 0;
			}
		}
	}
	return 1;
}

int searchRoutePrim(int s, int lp, bool isSetUp)
{
	int i, j, k;
	bool unconnectedFlag;
	int a, b;
	int to, cost;

	priority_queue<Node, vector<Node>, greater<Node> > q;

	Node Nodes[NODE_NUM];
	Node doneNode;

	for(i=0; i<NODE_NUM; i++){
		for(j=0; j<NODE_NUM; j++){
			if(link[i][j]>LINK_NUM) continue;
			unconnectedFlag = 0;
			for(k=0; k<lp_size[lp]; k++){
				if(spec[s+k][link[i][j]] == 1 || path_back[link[i][j]][lp]){
					unconnectedFlag = 1;
				}
			}
			if(!unconnectedFlag){
				Nodes[i].edges_to.push_back(j);
				Nodes[i].edges_cost.push_back(1);

			}
		}
	}

	for(i=0; i<NODE_NUM; i++){
		Nodes[i].done = false;
		Nodes[i].cost = -1;
		Nodes[i].nodeNum = i;
		Nodes[i].from = INF;
	}

	Nodes[source[lp]].cost = 0;

	q.push(Nodes[source[lp]]);

	while(!q.empty() && Nodes[dest[lp]].from == INF){
		doneNode = q.top();
		q.pop();
		if(doneNode.done) continue;
	  	doneNode.done = true;
	  	for(i = 0; i < doneNode.edges_to.size(); i++){
	    	to = doneNode.edges_to[i];
	    	cost = doneNode.cost + doneNode.edges_cost[i];
			    if(Nodes[to].cost < 0 || cost < Nodes[to].cost){
			    	Nodes[to].cost = cost;
					Nodes[to].from = doneNode.nodeNum;
					q.push(Nodes[to]);
				}
		}
	}

	int dest_node = dest[lp];
	int from_node = Nodes[dest_node].from;
	int hop_counter = 0;
	while(from_node < NODE_NUM && dest_node != from_node) {
	    dest_node = from_node;
	    from_node = Nodes[dest_node].from;
	    hop_counter++;
	}

	if (isSetUp)
	{
		limit_hop_prim[lp] = hop_counter + ADDITIONAL_HOP;
	} 
	else 
	{
		if (hop_counter > limit_hop_prim[lp])
		{
			return 0;
		}
	}

	a = dest[lp];
	b = Nodes[a].from;

	if (b < NODE_NUM) {
		for (j = 0; j < NODE_NUM; j++) {
			for (k = 0; k < NODE_NUM; k++) {
				if(link[j][k]>LINK_NUM)continue;
				path_prim[link[j][k]][lp] = 0;
			}
		}
		while (b < NODE_NUM && a != b) {
			path_prim[link[b][a]][lp] = 1;
			a = b;
			b = Nodes[a].from;
		}
		return 1;
	}else{
		return 0;
	}
}

int searchRouteBack(int s, int lp, bool isSetUp)
{
	int i, j, k;
	bool unconnectedFlag;
	int a, b;
	int to, cost;

	priority_queue<Node, vector<Node>, greater<Node> > q;

	Node Nodes[NODE_NUM];
	Node doneNode;

	for(i=0; i<NODE_NUM; i++){
		for(j=0; j<NODE_NUM; j++){
			if(link[i][j]>LINK_NUM) continue;
			unconnectedFlag = 0;
			for(k=0; k<lp_size[lp]; k++){
				if(spec[s+k][link[i][j]] == 1 || path_prim[link[i][j]][lp]){
					unconnectedFlag = 1;
				}
			}
			if(!unconnectedFlag){
				Nodes[i].edges_to.push_back(j);
				Nodes[i].edges_cost.push_back(1);
			}
		}
	}

	for(i= 0; i < NODE_NUM ; i++){
	    for(j=0;j < Nodes[i].edges_to.size() ; j++){
	    }
 	}

	for(i=0; i<NODE_NUM; i++){
		Nodes[i].done = false;
		Nodes[i].cost = -1;
		Nodes[i].nodeNum = i;
		Nodes[i].from = INF;
	}

	Nodes[source[lp]].cost = 0;

	q.push(Nodes[source[lp]]);

	while(!q.empty() && Nodes[dest[lp]].from == INF){
		doneNode = q.top();
	 	q.pop();
		if(doneNode.done) continue;
	  	doneNode.done = true;
	  	for(i = 0; i < doneNode.edges_to.size(); i++){
	    	to = doneNode.edges_to[i];
	    	cost = doneNode.cost + doneNode.edges_cost[i];
	    	if(Nodes[to].cost < 0 || cost < Nodes[to].cost){
	    		Nodes[to].cost = cost;
				Nodes[to].from = doneNode.nodeNum;
				q.push(Nodes[to]);
			}
		}
	}

	int dest_node = dest[lp];
	int from_node = Nodes[dest_node].from;
	int hop_counter = 0;
	while(from_node < NODE_NUM && dest_node != from_node) {
	    dest_node = from_node;
	    from_node = Nodes[dest_node].from;
	    hop_counter++;
	}

	if (isSetUp)
	{
		limit_hop_back[lp] = hop_counter + ADDITIONAL_HOP;
	} 
	else 
	{
		if (hop_counter > limit_hop_back[lp])
		{
			return 0;
		}
	}

	a = dest[lp];
	b = Nodes[a].from;

	if (b < NODE_NUM && a != b) {
		for (j = 0; j < NODE_NUM; j++) {
			for (k = 0; k < NODE_NUM; k++) {
				if(link[j][k]>LINK_NUM)continue;
				path_back[link[j][k]][lp] = 0;
			}
		}
		while (b < NODE_NUM && a != b) {
			path_back[link[b][a]][lp] = 1;
			a = b;
			b = Nodes[a].from;
		}
		return 1;
	}else{
		return 0;
	}
}

void delLp(int lp, int p)
{
	int index = ind_prim[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];


	if(a){
		if(p==0 || p==1){
			for(i=0; i<b; i++){
				for(j=0;j<LINK_NUM;j++){
					if(spec[index+i][j] == 0 && path_prim[j][lp] == 1) throw "プライマリパス消去エラー";
					spec[index+i][j] = path_prim[j][lp] ^	spec[index+i][j];
				}
			}
		}

		if(p==0 || p==2){
			index = ind_back[lp];
			if(index == INF) return;
			for(i=0; i<b; i++){
				for(j=0;j<LINK_NUM;j++){
					if(spec[index+i][j] == 0 && path_back[j][lp] == 1) throw "バックアップパス消去エラー";
					spec[index+i][j] = path_back[j][lp] ^ spec[index+i][j];
				}
			}
		}
	}
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