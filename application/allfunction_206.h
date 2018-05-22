/* This program is a simulation of hitless retuning */

/* If (algCall==0 || == 2) no retuning option */
/* When using middle fit, retune() after adding lp*/
/* Add retuning time for 1+1 */
/* Stop retuning if new request */
/*include sorted list*/
// Adding both primary and backup to the same list with initial state given
/* Use T_temp, put high value if not needed */
/* include sorted list*/
/* Including backup first */
/* Parallel realloc used for conventional, Check retuning function*/

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
//名前を簡潔に指定するための記述

// Declaring functions
//{
	int initialize(void);
	int reInitialize(void);
	int initializeEvent(void);
	int readInput(void);
	int readDemands(void);
	int firstFit(int);
	int lastFit(int);
	int checkFirstFit(int);
	int checkFirstFitRerouting(int);
	int checkLastFit(int);
	int retuneDown();
	int retuneUp();
	int retuneDown_0();
	int retuneDownNonRerouting_0();
	int retuneDownRerouting_0();
	int retuneUp_0();
	void retune();
	void retune_0();
	int removeLP(int);
	int removeLP_0(int);
	int asign(int, int);
	int asignRerouting(int,int);
	int printSpec(void);
	int minIndAllocLP(int);
	int minFragAllocLP(int);
	double checkFrag(int, int);
	int genDemands(void);
	int midFitAllocLP(int);

	void delFromList(int, int);
	void addToList(int, int, int);			//With list, lp index, sorting value

	int checkMoveUp(int);
	int checkMoveDown(int);
	void fixMiddle();

	int checkExactFit(int);
	int checkExactBp(int);
	int checkExactFitRerouting(int);
	int checkExactBpRerouting(int);

	int getPrimRoot(int, int);
	int getBackRoot(int, int);

	void deleteLP(int, int);
	void deleteLPRerouting(int, int);

	int firstFit1(int);
	int minIndAllocLP1(int);
	int minFragAllocLP1(int);
	int midFitAllocLP1(int);

	int rpAllocLP(int);

	int firstlastLP(int);

	int checkEOL(void);
	void updateSpecStat(void);

	int hopretune(int);

	int hopretune(int);
	int removeLP1_1(int, int);
	int firstFit1_1(int, int);
	int checkFirstBp(int);
	int checkFirstBpRerouting(int);
	int asignBp(int, int);
	int asignBpRerouting(int, int);
	int retuneBp(void);
	int removeBp(int);
	int reAllocBp(int);
	int readResult();
	int readResultPy();
	int readResultReroutingPy();
	void statDefrag(void);
	void statDefragPy(void);
	void statDefragReroutingPy(void);
	int statAlgo(void);
	int statReroutingAlgo(void);
	int readInput(void);
	int writeOutput(void);
	int writeOutputPy(void);
	int writeOutputReroutingPy(void);
	void addToList2(int, int, int);

	void delFromList2(int, int, int);
	double sum(double, int);
	double ave(double, int);
	double var(double, int);
	double standard(double, int);

//}

// Declaring global vars
//{
	int lp_ind=0, lp_size[M], source[M], dest[M], spec_ind[M];    // Each LP has an arrival index, a couple s/d,
															// and a given size. It is allocated to a spectrum index
	bool spec[S][L], path[N][N][L];				// Spectrum is a slots*links matrix, path(s,d,link) 1 if sd use l
	bool linked_path[N][N][N][N], linked_bp[N][N][N][N], linked_crosspath[N][N][N][N];				// Paths sharing at least a link

	int low_ind, high_ind;						// Lowest used index for lastfit and highest for firstfit
	int blocked;					// Counting blocked request //ブロッキングを起こしたリクエストの数
	int isactive[M];							// Tracking if lp is allocated and where (0-off, 1-firstf, 2-lastf)
	int t_req[M], t_hold[M], t_exp[M];			// Tracking request arrival time and lp expiration time
	double t_req_event[M], t_hold_event[M], t_exp_event[M];

	int hops[N][N], bhops[N][N];
	int link[N][N];

	int max_hold, last_lp;

	int midlp = M;

	unsigned seed1 = 123;      //time(NULL);						// Initializing generator
	unsigned seed2 =  125;

	struct lpNode{
		int x;							//lp index
		int y;    				  //lp length or size or both, for sorting purpose
		int z;							//Initial lp state, 1 prime, 0 for backups
		lpNode *next;
	};

	lpNode *activeList = NULL;
	lpNode *midFitList = NULL;		// This won't change, or we would lose the list in memory
	lpNode *backupList = NULL;
	lpNode *mixtList = NULL;
	lpNode *realList = NULL;
	lpNode *tempList = NULL;

	struct Node {
		// このノードから伸びるエッジの情報
		vector<int> edges_to;    // 各エッジの接続先のノード番号
	 	vector<int> edges_cost;  // 各エッジのコスト

		// ダイクストラ法のためのデータ
		bool done;  // 確定ノードか否か
		int cost;   // このノードへの現時点で判明している最小コスト
		int nodeNum;
		int from;

		//比較演算子のオーバーロード
		bool operator> (const Node &node) const {
    	return (cost > node.cost);
		}
	};

	bool path_rr[L][M], bp_rr[L][M];//リルーティングのための変数

	int part[N][N];

	int algoCall;
	int retOp;
	int eol_count;

	int t_temp =0;
	int last_ret= 0;
	int temp_max = T_temp;

	int A= A1;

	double t;

	int bp_ind[M];
	bool bp[N][N][L];

	int togOp;
	int realOp;
	int rerouteOp; //the number of rerouting operations

	bool lpState[M];
	bool bpState[M];

	int last_blocked = 0;   // Put to 1 if last allocation try was denied

	struct Event {			// path coming, path disrupution, defrag start
		
		double time;   // the time of event
		int lpNum;  // the lp number
		int type;  // 0:coming, 1:disruption

		//比較演算子のオーバーロード
		bool operator> (const Event &event) const {
    	return (time > event.time);
		}
	};

	priority_queue<Event, vector<Event>, greater<Event> > eventQueue; // 優先度付き待ち行列

	Event startEvent[M];
	Event endEvent[M];
	Event nowEvent;
	Event nextEvent;
	int defragCount;
	vector<Event> defragEvent;

//}

int retuneBp()
{
//	if(t == 700) printSpec();
	if(!algoCall) return 0;

	if(algoCall==1){
	//	statDefrag();
	//	cout << "checkpoint6 "<< endl;
		statAlgo();
		retuneDownNonRerouting_0();
	//	cout << "checkpoint7 "<< endl;
		return 0;
	}

	if(algoCall==2){
	//	statDefrag();
	//	cout << "checkpoint6 "<< endl;
		statReroutingAlgo();
		retuneDownRerouting_0();
	//	cout << "checkpoint7 "<< endl;
		return 0;
	}

	if(algoCall==3){
	//	statDefrag();
	//	cout << "checkpoint6 "<< endl;
		statDefragPy();
	//	cout << "checkpoint7 "<< endl;
		return 0;
	}

	if(algoCall==4){
	//	statDefrag();
	//	cout << "checkpoint6 "<< endl;
		statDefragReroutingPy();
	//	cout << "checkpoint7 "<< endl;
		return 0;
	}

	if(algoCall == 5){
		removeBp(0);
//		cout << " CHECKING 1" << endl;
		hopretune(0);
		reAllocBp(1);
//		cout << " CHECKING 2" << endl;
//		cout << " CHECKING 3" << endl;
	}

	if(algoCall == 6){
		removeBp(0);
		reAllocBp(1);
	}

	if(algoCall == 7) {
		removeBp(0);
		retuneDown_0();
		reAllocBp(1);
	}
	return 0;
}

// 合計の計算
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


void statDefrag(){
	writeOutput();
//	printSpec();
	system("ampl smpe_f.run > runresult.txt");
//	system("ampl smpe_fhm.run > runresult.txt");
	readResult();
//	printSpec();
}

void statDefragPy(){

	writeOutputPy();
	// printSpec();
	system("python ssr_lno.py");
//	system("ampl smpe_fhm.run > runresult.txt");
	readResultPy();
//	printSpec();
}

void statDefragReroutingPy(){

	writeOutputReroutingPy();
	// printSpec();
	system("python ssrr_lno.py");
//	system("ampl smpe_fhm.run > runresult.txt");
	readResultReroutingPy();
//	printSpec();
}

int readResult()
{
	int i,j;
	int a,b, lp;
	ifstream fin, fin1;

	for(i=0;i<S;i++){
		for(j=0;j<L;j++)  spec[i][j]= 0;
	}

	fin.open ("smpe_fhm_result.txt");
		if (!fin){
			cout <<"Cannot read fhm result file" << endl;
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

int readResultPy()
{
	int i,j;
	int a,b, lp;
	ifstream fin, fin1;

	for(i=0;i<S;i++){
		for(j=0;j<L;j++)  spec[i][j]= 0;
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

	for(i=0;i<S;i++){
		for(j=0;j<L;j++)  spec[i][j]= 0;
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
			for(i=0;i<L;i++){
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


int hopretune(int index)   				// Index is the starting point
{
	int a, b, lp, lp2;
	int ret_time = 0;
	int realMov = 1;
	int s1, d1, s2, d2;

	lpNode *cur;
	lpNode *cur1;
	lpNode *cur2;

	tempList = NULL;
	cur = mixtList;//lpnodeのカーソル
	while ( cur != NULL ) {
		addToList2(2, cur->x, cur->z);
	//	cout << "cur x, z= " << cur->x << ", " << cur->z << endl;
		cur= cur->next;
	}
	cur = tempList;
//	cout << "Before"  << endl ;
//	printSpec();

	realList = NULL;
	cur1 = NULL;
//	cout << "cur = mixtList " << endl<< endl;

	while ( cur != NULL ) {						// Checking all active LPs
		realList = NULL;
		cur1 = cur;
		cur = cur->next;
//		cout << "cur1 = cur " << endl<< endl;

		while ( cur1 != NULL ) {	  // List for parallel realloc.
			lp = cur1->x;
		//	int st = cur1->z;
			cur1 = cur1->next;

			cur2 = realList;
//			cout << "cur2 = realList 1 " << endl<< endl;

			if(cur2 == NULL){
				addToList2(1, lp, 1);
				if(cur && cur->x == lp) cur = cur->next;
				delFromList2(4, lp, 1);
			}else{
				int conf = 0;
				while(cur2 != NULL){
					s2 = source[cur2->x];
					d2= dest[cur2->x];
					cur2 = cur2->next;
					s1 = source[lp];
					d1= dest[lp];
					if(linked_path[s1][d1][s2][d2] || linked_bp[s1][d1][s2][d2] || linked_crosspath[s1][d1][s2][d2]){
						conf =1;
						break;
					}
				}
				if(conf==0){
					addToList2(1, lp, 1);
					if(cur && cur->x == lp) cur = cur->next;
					delFromList2(4, lp, 1);
				}
			}
		}

		cur2 = realList;
	//	cout << "cur2 = realList 2" << endl<< endl;
		int realcheck = realOp;
		while ( cur2 != NULL ){				// Parallel realloc.
			lp = cur2->x;
		//	b =  cur2->z;
			cur2 = cur2->next;
			if(spec_ind[lp] > index){
				deleteLP(lp, 1);        // Remove LP from spec to avoid self-conflict
			//	printSpec();
				a = checkFirstFit(lp);				// Return spec_ind[lp] at worst
				if(a==S || a > spec_ind[lp]) a = checkFirstFit(lp);
	//			cout << "In while-If"  << endl ;
				if(a != spec_ind[lp]){
					togOp++;
					realOp++;
					realMov++;
					lpState[lp]= 0;
					bpState[lp]= 1;
				} // retOp++;
				if(a< INF) asign(lp, a);
			//	if(a>=INF) blocked++;
	//			cout << "Out while-If"  << endl ;
			//	ret_time++ ;
			//	if(t+ret_time >= t_req[last_lp+1] || ret_time >= temp_max) return 0;
			}
		}
		cur2 = realList;
		while ( cur2 != NULL ){
			lp = cur2->x;
			cur2 = cur2->next;
			delFromList(3, lp);
		}
		if(realcheck != realOp){
			ret_time++ ;
			if(t+ret_time >= t_req[last_lp+1] || ret_time >= temp_max) return 0;
		}
	}
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
//	cout << "cur = realList " << endl<< endl;
	// while ( cur != NULL ) {
		// cout << "cur = " << cur->x <<", "<<cur->z << endl;
		// cur= cur->next;
	// }

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
					if(a==S || a > spec_ind[lp]) a = checkFirstFit(lp);
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
					if(a==S || a > spec_ind[lp]) a = checkFirstBp(lp);
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
				ret_time += 1/double(K);					//increment defrag time
				nextEvent = eventQueue.top(); 	// next event
				// cout << "nextEvent.type" << nextEvent.type << endl;
				if(nextEvent.type == 0 && (t+ret_time >= nextEvent.time || ret_time >= temp_max)){
					// t += ret_time;
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
//	cout << "cur = realList " << endl<< endl;
	// while ( cur != NULL ) {
		// cout << "cur = " << cur->x <<", "<<cur->z << endl;
		// cur= cur->next;
	// }

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
								for(i=0;i<L;i++){
									if(path_rr[i][cur2->x] && path_rr[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<L;i++){
									if(path_rr[i][cur2->x] && bp_rr[i][lp])conf =1;
								}

							}
						}else{//バックアップ
							if(lpState[lp]){
								for(i=0;i<L;i++){
									if(bp_rr[i][cur2->x] && path_rr[i][lp]) conf =1;
								}
							}else{
								for(i=0;i<L;i++){
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
					deleteLPRerouting(lp, 1);// Remove LP from spec to avoid self-conflict
					//プライマリパスのみ消去
					a = checkExactFitRerouting(lp);
					// std::cout << " after checkExactFitRerouting a = " << a << '\n';
					if(a==S || a > spec_ind[lp]){
						a = checkFirstFitRerouting(lp);
						// std::cout << " after checkFirstFitRerouting a = " << a << '\n';
					}
					//もしスロット番号が小さくならないようであれば
				//	if(a > spec_ind[lp]) a = spec_ind[lp];
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
					// std::cout << "バックアップ, lp = " << lp << ", bp_ind[lp] = " << bp_ind[lp] << ", source[lp] = " << source[lp] << ", dest[lp] =" << dest[lp] << ", lp_size[lp] =" <<lp_size[lp]<<'\n';
					deleteLPRerouting(lp, 2);            // Remove LP from spec to avoid self-conflict
					//バックアップパスのみ消去
					a = checkExactBpRerouting(lp);
					// std::cout << " after checkExactFitRerouting a = " << a << '\n';
					if(a==S || a > bp_ind[lp]){
						a = checkFirstBpRerouting(lp);
						// std::cout << " after checkFirstFitRerouting a = " << a << '\n';
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
				ret_time += 1/double(K);					//increment defrag time
				nextEvent = eventQueue.top(); 	// next event
				// cout << "nextEvent.type" << nextEvent.type << endl;
				if(nextEvent.type == 0 && (t+ret_time >= nextEvent.time || ret_time >= temp_max)){
					// t += ret_time;
					return 0;	
				}
				//新しいパスがきているか, 100000ステップを超えたならばデフラグ終了
			}

		}//tempListがなくなるまで
	}//realMovが0以外なら続くのwhile文 帯域移動がなければデフラグ終了
	return 0;
}

int retuneDown_0()
{
	int s, d;
	int index1, index2;
	int i,j,p;
  	bool eol=0;
	int k=0;

	int ret_time = 0 ;
	int mov_time = 0 ;
//	int ret_time = t ;

	index1 = 1;

	while(index1 < S-1 && k++ < INF){			//Sweep the spectrum
		index2 = index1-1;
		lpNode *cur = activeList; 
		mov_time = 0 ;
		while ( cur != NULL ) {						// Checking all active LPs
			i = cur->x;
			cur = cur->next;
			if(isactive[i] == 1 && spec_ind[i] == index1){
			//	cout << "1: i " << i << ", Spec ind "<< spec_ind[i] << " index1 " << index1 <<  ",  index2 " << index2 << endl;
				s = source[i];
				d = dest[i];
				eol = 0;
				while(index2 >= 0 && !eol){		// Check to what extend it can be retuned

					for(j=0;j<L;j++){							//Check path availibility
						if(spec[index2][j] && path[s][d][j]){
						//	cout << "2: i " << i << ", Spec ind "<< spec_ind[i] << " index1 " << index1 <<  ",  index2 " << index2 << endl;
							eol = 1;
							break;
						}
					}
					if(eol) break;
					index2--;
				}
			}

		//	cout << "3: i " << i << ", Spec ind "<< spec_ind[i] << " index1 " << index1 <<  ",  index2 " << index2 << endl;

			if(isactive[i] &&(index2 + 1) < index1 && spec_ind[i] == index1){			// Retune retunables
				for(p=0; p<lp_size[i]; p++){						// an xor with its path will remove it
					for(j=0;j<L;j++){
						spec[index1+p][j] = path[s][d][j] ^	spec[index1+p][j];
					}
				}

				asign(i, (index2 + 1));

				if(index1 - index2 > mov_time) mov_time = index1 - index2;

				index2= index1 - 1;
				retOp++;

		//		cout << "LP index " << i << " retuned down from index " << index1 <<", New index " << spec_ind[i] << endl;	//Just for checking
		//		printSpec();
			}
		}
		ret_time += mov_time;
		if(t+ret_time >= t_req[last_lp+1]) return 0; //if nextpath comming
		if(last_blocked && ret_time >= temp_max) return 0; 
		index1++;
	}
	return 0;
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

	while(index1 < S-1 && k++ < INF){			//Sweep the spectrum
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
				if(a==S || a > spec_ind[i]) a = checkFirstFit(i);
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
		ret_time += double(mov_time)/double(K);
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

	while(index1 < S-1 && k++ < INF){			//Sweep the spectrum
		lpNode *cur = activeList; 
		mov_time = 0 ;
		while ( cur != NULL ) {						// Checking all active LPs
			i = cur->x;
			cur = cur->next;
			if(isactive[i] == 1 && spec_ind[i] == index1){
				deleteLPRerouting(i, 1);// Remove LP from spec to avoid self-conflict
				//プライマリパスのみ消去
				a = checkExactFitRerouting(i);
				// std::cout << " after checkExactFitRerouting a = " << a << '\n';
				if(a==S || a > spec_ind[i]){
					a = checkFirstFitRerouting(i);
					// std::cout << " after checkFirstFitRerouting a = " << a << '\n';
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
				deleteLPRerouting(i, 2);            // Remove LP from spec to avoid self-conflict
				//バックアップパスのみ消去
				a = checkExactBpRerouting(i);
				// std::cout << " after checkExactFitRerouting a = " << a << '\n';
				if(a==S || a > bp_ind[i]){
					a = checkFirstBpRerouting(i);
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
		ret_time += double(mov_time)/double(K);
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
				for(j=0;j<L;j++){//全てのリンクについて
					spec[index+i][j] = path_rr[j][lp] ^ spec[index+i][j];
					//全てのリンクの占有する帯域スロット番号について
					//プライマリパスが通っているところはビット反転させる
				}
			}

			index = bp_ind[lp];
			if(index < INF){
				for(i=0; i<b; i++){									// an xor with its path will remove it
					for(j=0;j<L;j++){
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
				for(j=0;j<L;j++){//全てのリンクについて
					spec[index+i][j] = path[s][d][j] ^ spec[index+i][j];
					//全てのリンクの占有する帯域スロット番号について
					//プライマリパスが通っているところはビット反転させる
				}
			}

			index = bp_ind[lp];
			if(index < INF){
				for(i=0; i<b; i++){									// an xor with its path will remove it
					for(j=0;j<L;j++){
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

	while(index <= (S-b) && !asigned)   		  //Checking all spectrum range
	{
		if(!b) break;

		for(i=0;i<L;i++){							//Check path availibility
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
		if(i==L && !nofit){
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

	while(index <= (S-b) && !asigned)   		  //Checking all spectrum range
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


//	if(t >= 1300) cout << "1: high_ind for lp, ind=" << lp <<", "<< index << endl;

	bp_ind[lp] = index;

//	if(t >= 1300) cout << "2: high_ind at t=" << t << endl;

	for(j=0;j<b;j++){
		for(p=0;p<L;p++){
			if(spec[index+j][p] == 1 && bp[s][d][p] ==1) throw "バックアップパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || bp[s][d][p];
		}
	}
	return 0;
}

int asignBpRerouting(int lp, int index)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];


//	if(t >= 1300) cout << "1: high_ind for lp, ind=" << lp <<", "<< index << endl;

	bp_ind[lp] = index;

//	if(t >= 1300) cout << "2: high_ind at t=" << t << endl;
	// for(i = 0; i < L ; i ++){
	// 	cout << "bp_rr[" << i << "][" << lp << "] = " << bp_rr[i][lp] << endl;
	// }
	for(j=0;j<b;j++){
		for(p=0;p<L;p++){
			if(spec[index+j][p] == 1 && bp_rr[p][lp] ==1) throw "バックアップパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || bp_rr[p][lp];
		}
	}

	return 0;
}

int removeBp(int a)   // a=0 for backups and 1 for primaries
{
	int lp;

	lpNode *cur = activeList;//primaryのlist

	while ( cur != NULL ){						// Checking all active LPs
		lp = cur->x;
		cur = cur->next;
		if(!a) deleteLP(lp, 2); // Remove backup without desactivating primary
		if(a) deleteLP(lp, 1);  // Remove primary without desactivating backup
	}
	return 0;
}

int reAllocBp(int p)
{
	int lp, a, index=0;

	if(p){
		lpNode *cur = activeList;
		while ( cur != NULL ) {						// Checking all active LPs
			lp = cur->x;
			cur = cur->next;
		//	if(t >= 600 && lp==33) cout << "1: high_ind = " << bp_ind[lp] << endl;
			a = checkFirstBp(lp);
		//	if(t >= 600 && lp==33){
		//		cout << "2: Realloc at lp" <<lp <<" at" << a << endl;
		//		printSpec();
		//	}
			if(bp_ind[lp] != a){
				realOp++;
			//	if(bpState[lp]) togOp++;
			}

			if(a<INF) asignBp(lp, a);
			bp_ind[lp] = a;
		}
	}

	if(!p){
		while(index < S){			//Sweep the spectrum
			lpNode *cur = activeList;
			while ( cur != NULL ) {						// Checking all active LPs
				lp = cur->x;
				cur = cur->next;
				if(bp_ind[lp] == index){
					// if(t >= 700){
						// cout << "1: high_ind for lp, ind= " << lp <<", "<< index << endl;
						// if(lp == 181) printSpec();
					// }
					a = checkFirstBp(lp);
				//	if(t >= 700) cout << "2:  high_ind for lp, ind= " << lp <<", "<< a << endl;
					if(a < INF) asignBp(lp, a);
					bp_ind[lp] = a;
				}
			}
			index++;
		//	if(t >= 1300) cout << "3: high_ind at t=" << t << endl;
		}
	}
	return 0;
}

int initialize(void)						  // Set everything to zero
{
	int i,j, p, k, l;

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			for(p=0;p<L;p++) path[i][j][p] = 0;
			hops[i][j]=0;
			part[i][j]=0;
			link[i][j]= INF;
			for (k=0;k<N;k++){
				for (l=0;l<N;l++){
					linked_path[i][j][k][l] = 0;
					linked_bp[i][j][k][l] = 0;
					linked_crosspath[i][j][k][l] = 0;
				}
			}
		}
	}

	for(i=0; i<L; i++){
		for(j=0;j<M;j++){
			path_rr[i][j] = 0;
			bp_rr[i][j] = 0;
		}
	}

	reInitialize();
	return 0;
}

int reInitialize(void)						// Set everything to zero save routings and demands
{
	int i,j;

	low_ind = S-1; high_ind=0; lp_ind=0;
	blocked=0;
	retOp = 0;
	eol_count = 0;
	togOp =0;
	realOp = 0;
	rerouteOp = 0;

	for(i=0;i<M;i++){
		spec_ind[i]=0; isactive[i]=0;
		lpState[i]=1; bpState[i]=0;
	}

	//make priority empty
	while(!eventQueue.empty()){
		eventQueue.pop();
	}

	for(i=0;i<S;i++){
		for(j=0;j<L;j++)  spec[i][j]= 0;
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

	for(i=0;i<M;i++){
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

int readDemands()
{
	int l=0;
	char tmp1='r';
	int i;

	ifstream fin;
					//Count number of demands from file

	fin.open ("input_demands.txt");        //Input demands from file
		if (!fin)
		{
			cout <<"Cannot open demands file" << endl;
			return 1;
		}

		fin.ignore(INT_MAX,'=');
		fin.ignore(INT_MAX,':');
		fin >> tmp1;					//tmp1 is a char

		while(tmp1 !=';' && l < M){
			l++;
			fin.ignore(INT_MAX,':');
			fin >> tmp1;
		}
	fin.close();

		cout << "Nbr of demands := " << l << endl;

		//l1 = l;
		//if(l==99) return 0;

					// Read demands

	fin.open ("input_demands.txt");
		if (!fin){
			cout <<"Cannot open demands file" << endl;
			return 1;
		}

		fin.ignore(INT_MAX,'=');

		for (i=0;i<l;i++){
			fin.ignore(INT_MAX,':');
			fin >> source[lp_ind] >> dest[lp_ind] >> lp_size[lp_ind] >> t_req[lp_ind] >> t_hold[lp_ind];
			t_exp[lp_ind] = t_req[lp_ind] + t_hold[lp_ind];
			if(source[lp_ind] != dest[lp_ind])
				lp_ind++;
		}
	fin.close();
	return 0;
}

int genDemands()
{
	int i,j,k;
	int arr_int;
	double arr_int_event;
	double temp1;

	double mu = 1/double(H);
	// for exponentially-distributed duration with mean H=1/mu
	//指数関数分布の継続時間 0.1
	double inter_arr = double(H)/A;
	// Inter-arrival time for poisson distribution A=H/inter_arr
	//ポアソン分布の到着間隔 10/74

	srand (seed2);					//time(NULL)2種類ある

	default_random_engine generator (seed1);
	//generatorという擬似乱数生成エンジンを定義 種はseed1
	// poisson_distribution<int> next_arr(K*inter_arr);
	// exponential_distribution<double> next_arr(1/(K*inter_arr));
	exponential_distribution<double> next_arr(1/inter_arr);
	// 100* to clock the steps to a 10ms, does not change A,
	//平均値K*inter_arr(100*10/74)でnext_arrというポアソン分布が得られる
	exponential_distribution<double> hold_time(mu);
	// A = 100*H/100*inter_arr
	// ρ = λ × H
	//   = 100 * 10 / 74 * 10
	//	 = 135
	// 平均値mu(0.1)でhold_timeという指数関数分布が得られる
	uniform_int_distribution<int> traff_dist(1, req_Max);
	//1からreq_Maxまでの整数がtraff_distで等確率で得られる

	ofstream ofs1;
    ofs1.open("./../result/input_demands1.txt");
	if(!ofs1){
		cout<< "Cannot open input_demands1 file"<<endl;
		return 1;
	}

	ofstream ofs2;
    ofs2.open("./../result/input_demands2.txt");
	if(!ofs2){
		cout<< "Cannot open input_demands2 file"<<endl;
		return 1;
	}

	ofs1 << "mu and inter:=" << mu <<", "<<  inter_arr<< endl;
	ofs1 << "LP number, s, d, size, Arrival time, holding time:=" << endl;
	ofs2 << "mu and inter:=" << mu <<", "<<  inter_arr << endl;
	ofs2 << "LP number, s, d, size, Arrival time, holding time:=" << endl;

	max_hold = 0;
	t_req[0]= 0;
	t_req_event[0] = 0;

	for (int i=0; i<M; ++i){ //リクエストの数だけ
		temp1 = hold_time(generator); //指数関数分布で継続時間を乱数として得る
		t_hold[i] = int(temp1)+1;		//100をかけて1を足したものが継続時間となる
		t_hold_event[i] = temp1;
		if(max_hold < t_hold[i]){
			max_hold = t_hold[i];//継続時間の最大値を超えていれば最大値とする
		}
		arr_int_event = next_arr(generator);
		arr_int = int(arr_int_event)+1; //ポアソン分布で到着間隔を得る
		lp_size[i] = traff_dist(generator); //等確率で占有帯域スロット数を得る
		source[i] = rand() %N;//発ノードをランダムで得る
		dest[i] = rand() %N;//着ノードをランダムに得る
			while(source[i] == dest[i]) dest[i] = rand() %N;//発ノードと着ノードが同じにならないようにする
		t_exp[i] = t_req[i] + t_hold[i]; //継続時間と到着時刻を合わせて切断時刻を計算
		t_exp_event[i] = t_req_event[i] + t_hold_event[i]; //event driven
		ofs1 << i  << ": "<< source[i]  <<" " << dest[i]  <<" " << lp_size[i]  <<" " << t_req[i] <<" " << t_hold[i] << endl;
		ofs2 << i  << ": "<< source[i]  <<" " << dest[i]  <<" " << lp_size[i]  <<" " << t_req_event[i] <<" " << t_hold_event[i] << endl;
		t_req[i+1]= t_req[i] + arr_int; //到着時刻と到着間隔を合わせて次の到着時刻を計算
		t_req_event[i+1] = t_req_event[i] + arr_int_event; //event driven
		//lp番号ごとのルートを持つ変数も定める
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				if(link[j][k] < L){
					path_rr[link[j][k]][i] = path[source[i]][dest[i]][link[j][k]];
					bp_rr[link[j][k]][i] = bp[source[i]][dest[i]][link[j][k]];
				}
			}
		}
	}

	ofs1 << ":;"<< endl << endl;
	ofs1.close();
	ofs2 << ":;"<< endl << endl;
	ofs2.close();
	return 0;
}

void retune_0()
{
	int i, j;
	int a;

	retuneDown_0();
	retuneUp_0();

	updateSpecStat();

	fixMiddle();
}

int retuneUp_0()
{
	int s, d;
	int index1, index2, a;
	int i,j,p;
  	bool eol=0;
	int k=0;

	index1 = S-2;

	while(index1 >= 0 && k++ < INF){			//Sweep the spectrum
		index2 = index1;

		lpNode *cur = activeList;
		while ( cur != NULL ) {						// Checking all active LPs
			i = cur->x;
			cur = cur->next;
			a = lp_size[i];
	//		if(isactive[i]&& spec_ind[i]+a == index1)		// && spec_ind[i]+a == index1
	//			cout << "isactive " << isactive[i] <<", Spec ind " << spec_ind[i] << ", index1 " << index1<< endl;
	//		cout << "1: i " << i << " of size "<< a << ", Spec ind "<< spec_ind[i]+a << " index1 " << index1 <<  ",  index2 " << index2 << endl;
			if(isactive[i]==2 && spec_ind[i] + a == index1){
				s = source[i];
				d = dest[i];
				eol = 0;

				while(index2 < S && !eol){		// Check to what extend it can be retuned
					for(j=0;j<L;j++){							//Check path availibility
						if(spec[index2][j] && path[s][d][j]){

	//		cout << "2: i " << i << ", Spec ind "<< spec_ind[i]+a << " index1 " << index1 <<  ",  index2 " << index2 << endl;
							eol = 1;
							break;
						}
					}
					if(eol) break;
					index2++;
				}
			}

	//		cout << "i " << i << ", Spec ind " << spec_ind[i] << ", index1 = " << index1 <<  ",  index2 = " << index2 << endl;

			if(isactive[i] && index2 > index1 && spec_ind[i] + a == index1){						// Retune retunables
			//	removeLP(i);
	//		cout << "3: i " << i << ", Spec ind "<< spec_ind[i]+a << " index1 " << index1 <<  ",  index2 " << index2 << endl;
				for(p=0; p<lp_size[i]; p++){									// an xor with its path will remove it
					for(j=0;j<L;j++){
						spec[spec_ind[i]+p][j] = path[s][d][j] ^ spec[spec_ind[i]+p][j];
					}
				}
				int c= spec_ind[i];
				isactive[i] = 2;
				asign(i, (index2 - a));
				index2= index1;	 		//index2= index1 + 1
				retOp++;

	//			cout << "LP index " << i <<" of size " << lp_size[i]<< " retuned up from " << c <<" , New index " << spec_ind[i] << endl;	//Just for checking
			//	printSpec();
			}
		}
		index1--;
	}
	return 0;
}

void retune(){
	int i, j, a, b;

	if(low_ind - high_ind - 1 < 4*req_Max){
//		releaseMiddle();
		isactive[midlp] = 1;
		midlp = M;
	}

	if(t_temp <= 0){						//	Retune only after the previous one is done
		a = retuneDown();
//	if(t==T_end-3) cout << " 3_0 lp index " << spec_ind[9989] <<endl;
		b = retuneUp();
//	if(t==T_end-3) cout << " 3_1 lp index " << spec_ind[9989] <<endl;

		updateSpecStat();

		fixMiddle();

		if(a > b){	t_temp = a;
		}else t_temp = b;
	}
}

int retuneDown()
{
	int s, d;
	int index1, index2;
	int i,j,p;
  	bool eol=0;
	int k=0;

	index1 = 1;

//	if(t>=T_end-3) cout << " 3_0 lp index " << spec_ind[9989] <<endl;

	while(index1 < S-1 && k++ < INF){			//Sweep the spectrum
		index2 = index1-1;
		lpNode *cur = activeList;
		while ( cur != NULL ) {						// Checking all active LPs
			i = cur->x;
			cur = cur->next;
		//	if(i == 9989) cout << " 3_1 lp index " << spec_ind[9989] <<endl;
			if(isactive[i] == 1 && spec_ind[i] == index1){
				s = source[i];
				d = dest[i];
				eol = 0;
			//	while(index2 >= 0 && !eol){		// Check to what extend it can be retuned

					for(j=0;j<L;j++){							//Check path availibility
						if(spec[index2][j] && path[s][d][j]){
							eol = 1;
							break;
						}
					}
			//		if(eol) break;
			//		index2--;
			//	}
			}

			if(isactive[i]==1 && !eol && spec_ind[i] == index1){		// Retune retunables
				for(p=0; p<lp_size[i]; p++){									// an xor with its path will remove it
					for(j=0;j<L;j++){
						spec[index1+p][j] = path[s][d][j] ^	spec[index1+p][j];
					}
				}
			//	if(t>=T_end-3) cout << " 3 lp index " << spec_ind[9989] << ", index2 " << index2 <<endl;
				asign(i, (index2));
				retOp++;
				last_ret = i;
			//	if(t>=T_end-3) cout << " 4 lp index " << spec_ind[9989] <<endl;
				return index1 - index2;
				index2= index1 - 1;
			}
		}
		index1++;
	}
	return 0;
}

int retuneUp()
{
	int s, d;
	int index1, index2, a;
	int i,j,p;
  	bool eol=0;
	int k=0;

	index1 = S-2;

	while(index1 >= 0 && k++ < INF){			//Sweep the spectrum
		index2 = index1;

		lpNode *cur = activeList;
		while ( cur != NULL ) {						// Checking all active LPs
			i = cur->x;
			cur = cur->next;
			a = lp_size[i];

			if(isactive[i]==2 && spec_ind[i] + a == index1){
				s = source[i];
				d = dest[i];
				eol = 0;

			//	while(index2 < S && !eol){		// Check to what extend it can be retuned
					for(j=0;j<L;j++){							//Check path availibility
						if(spec[index2][j] && path[s][d][j]){
							eol = 1;
							break;
						}
					}
			//		if(eol) break;
			//		index2++;
			//	}
			}
			if(isactive[i]==2 && !eol && spec_ind[i] + a == index1){						// Retune retunables
				for(p=0; p<lp_size[i]; p++){									// an xor with its path will remove it
					for(j=0;j<L;j++){
						spec[spec_ind[i]+p][j] = path[s][d][j] ^ spec[spec_ind[i]+p][j];
					}
				}
				int c= spec_ind[i];
				isactive[i] = 2;
				asign(i, (index2+1 - a));
				retOp++;
				last_ret = i;

				return index2+1 - a - index1;
				index2= index1;
			}
		}
		index1--;
	}
	return 0;
}

int checkEOL()
{
	int i, j;
	int a;

	lpNode *cur;

	cur = activeList;
	while ( cur != NULL ) {			// Checking active LPs and their position
		i = cur->x;
		cur = cur->next;
		if(isactive[i] == 1){				// Updating high
			a = checkFirstFit(i);
			if(a < spec_ind[i]){
				eol_count++;
			//	return 1;
			}
		}
		if(isactive[i] == 2){				// Updating low
			a = checkLastFit(i);
			if(a != INF && a > spec_ind[i]){
				eol_count++;
			//	return 1;
			}
		}
	}
	return 0;
}

void updateSpecStat()
{
	int i, j;
	int new_high= 0, new_low= S-1;								// Updating high and low index
	int attached = 0;

	high_ind =0;
	low_ind =S-1;

	lpNode *cur;

	cur = activeList;
	while ( cur != NULL ) {			// Checking active LPs and their position
		i = cur->x;
		cur = cur->next;
		new_high = spec_ind[i]+lp_size[i] -1;
		new_low = spec_ind[i];
		if(isactive[i] == 1 && new_high >  high_ind){				// Updating high
			high_ind = new_high;
		}
		if(isactive[i] == 2 && new_low <  low_ind){				// Updating low
			low_ind = new_low;
		}
	}
}

int firstFit(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = 0;


	if(!b) return 0;

	while(index <= (S-b) && !asigned)   		  //Checking all spectrum range
	{
	//	cout << "LP checking index: " << lp << ", "<< index << endl;

		for(i=0;i<L;i++){							//Check path availibility
			if(spec[index][i] && path[s][d][i]){
				index++;
				break;
			}else{
				nofit = 0;
				for(j=0;j<b;j++){							//Checking if it fit (size)
					if(spec[index+j][i] && path[s][d][i]){
						index += j;
						nofit = 1;
						break;
					}
				}
			}
			if(nofit) break;
		}
		if(i==L && !nofit){
			isactive[lp]= 1;			// Side is defined before asign
			asign(lp, index);
			asigned= 1;
	//		cout << "LP index " << lp << " first-asigned to index" << spec_ind[lp] << endl;
		//	cout << "LP index is allocated #" << isactive[lp] << endl;
		}
	}
	if(!asigned){
		blocked ++;
		return 0;
	}else return 1;
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

	while(index <= (S-b) && !asigned) //Checking all spectrum range 全ての帯域スロットを調べる
	{
		if(!b) break;//bが0ならブレイク

	//	cout << "LP checking index: " << lp << ", "<< index << endl;

		for(i=0;i<L;i++){		//Check path availibility
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
		if(i==L && !nofit){//もし全てのリンクを調べ終え、リンクが使われていなければ
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

	while(index <= (S-b) && !asigned) //Checking all spectrum range 全ての帯域スロットを調べる
	{
		if(!b) break;//bが0ならブレイク
		//	cout << "LP checking index: " << lp << ", "<< index << endl;
		isGetRoot = getPrimRoot(index, lp);
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

int lastFit(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = S-b;

	while(index >= 0 && !asigned)   		  //Checking all spectrum range
	{
		if(!b) break;

	//	cout << "LP checking index: " << lp << ", "<< index << endl;
		for(i=0;i<L;i++){							//Check path availibility
	//		cout << "before break" << endl;
			if(spec[index][i] && path[s][d][i]){
				index--;
				break;
			}else{
				nofit = 0;
				for(j=0;j<b;j++){							//Checking if it fit (size)
					if(spec[index+j][i] && path[s][d][i]){
						index -= b-j; 		//Jump to the next fitable index
						nofit = 1;
						break;
					}
				}
			}
			if(nofit)	break;
	//		cout << "Index " << index << ", lp " << lp << endl;
		}

		if(i==L && !nofit){
			isactive[lp]= 2;
			asign(lp, index);
			asigned= 1;
	//		cout << "LP index " << lp << " last-asigned to index " << spec_ind[lp] << endl;
		//	cout << "LP index is allocated #" << isactive[lp] << endl;
		}
	}
	//cout << "LP index assigned if 1: " << asigned << endl << endl;
	if(!asigned){
		blocked ++;
		return 0;
	}else return 1;
}

int checkLastFit(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = S-b;

	while(index >= 0 && !asigned)   		  //Checking all spectrum range
	{
		if(!b) break;

	//	cout << "LP checking index: " << lp << ", "<< index << endl;
		for(i=0;i<L;i++){							//Check path availibility
	//		cout << "before break" << endl;
			if(spec[index][i] && path[s][d][i]){
				index--;
				break;
			}else{
				nofit = 0;
				for(j=0;j<b;j++){							//Checking if it fit (size)
					if(spec[index+j][i] && path[s][d][i]){
						index -= b-j; 		//Jump to the next fitable index
						nofit = 1;
						break;
					}
				}
			}
			if(nofit)	break;
	//		cout << "Index " << index << ", lp " << lp << endl;
		}

		if(i==L && !nofit){
			asigned= 1;
	//		cout << "LP index " << lp << " will be last-asigned to index " << index << endl;
			return index;
		//	cout << "LP index is allocated #" << isactive[lp] << endl;
		}
	}
	//cout << "LP index assigned if 1: " << asigned << endl << endl;
	if(!asigned)
		return INF;
	return 0;
}

int firstlastLP(int lp)
{
	int a, b, c;

	int s = source[lp], d= dest[lp];

	c = lp_size[lp];
	if(!c) return 0;

//	if(lp == 561)
	//	cout << "0: LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;

	if(lp%2){
		a= checkFirstFit(lp);				// Available index for first or last fit

	//	if(lp == 561)
		//	cout << "1: LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;

		if(a < low_ind){
			isactive[lp]= 1;
			asign(lp, a);
//			cout << "LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}else{
			b= checkLastFit(lp);
			if(b == INF){			// No available index to fit in
				blocked ++;
				return 0;
			}else{
				isactive[lp] = 2;
				asign(lp, b);
				return 1;
			}
		}
//		if(lp == 561)
//			cout << "1: LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
	}else{
		b= checkLastFit(lp);
//		a=b;

//		if(lp == 561)
//			cout << "2: LP index " << lp << " last-asigned_0 to index  " << b << endl;

		if(b > high_ind && b != INF){
			isactive[lp]= 2;
			asign(lp, b);
//			cout << "LP index " << lp << " last-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}else{
			a= checkFirstFit(lp);
			if(a == INF){			// No available index to fit in
				blocked ++;
				return 0;
			}else{
				isactive[lp] = 1;
				asign(lp, a);
				return 1;
			}
		}
//		if(lp == 561)
//			cout << "2: LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
	}
}

int minIndAllocLP(int lp)
{
	int a, b, c;

	c = lp_size[lp];
	if(!c) return 0;

	a= checkFirstFit(lp);				// Available index for first or last fit
	b= checkLastFit(lp);

	if(a == INF && b == INF){			// No available index to fit in
		blocked ++;
		return 0;
	}

	{	// min max index Allocation policy
		if(b > low_ind && a+c < high_ind){			// No added index in both side
			if(high_ind > S-low_ind){ 		// Lower side higher
				isactive[lp]= 2;
				asign(lp, b);
	//			cout << "LP index " << lp << " last-asigned_0 to index  " << spec_ind[lp] << endl;
				return 1;
			}else{							// Upper side higher
				isactive[lp]= 1;
				asign(lp, a);
	//			cout << "LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
				return 1;
			}
		}else{							// Added index in at least one side
			if(a + c -1 -high_ind < low_ind - b){ 		// Added index lower for first fit
				isactive[lp]= 1;
				asign(lp, a);
	//			cout << "LP index " << lp << " first-asigned_1 to index " << spec_ind[lp] << endl;
				return 1;
			}
			if(a + c -1 -high_ind > low_ind - b){ 		// Added index lower for last fit
				isactive[lp]= 2;
				asign(lp, b);
	//			cout << "LP index " << lp << " last-asigned_1 to index " << spec_ind[lp] << endl;
				return 1;
			}
			if(a + c -1 - high_ind == low_ind - b && high_ind +1 > S-low_ind){ 		// Added same for first and last fit, but last fit shorter
				isactive[lp]= 2;
				asign(lp, b);
	//			cout << "LP index " << lp << " last-asigned_2 to index " << spec_ind[lp] << endl;
				return 1;
			}else{										// first-fit shorter or the same size
				isactive[lp]= 1;
				asign(lp, a);
	//			cout << "LP index " << lp << " first-asigned_2 to index " << spec_ind[lp] << endl;
				return 1;
			}
	//	cout << "isactive a=" << isactive[lp] << endl;
		}
	}
}

int minFragAllocLP(int lp)
{
	int a, b, c;

	c = lp_size[lp];
	if(!c) return 0;

	a= checkFirstFit(lp);				// Available index for first or last fit
	b= checkLastFit(lp);

	if(a == INF && b == INF){			// No available index to fit in
		blocked ++;
		return 0;
	}

	{									// min max index Allocation policy
		double frag1 = checkFrag(lp, a);
		double frag2 = checkFrag(lp, b);
	//	cout << "LP index " << lp << ", a " << a << ",and frag1 " << frag1<< ", b " << b << ",and frag2 " << frag2 << endl;

		if(frag1 < frag2){				// First fit causes less fragmentation
			isactive[lp]= 1;
			if(a > low_ind) isactive[lp]= 2;
			asign(lp, a);
	//		cout << "LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}
		if(frag2 < frag1){				// Last fit causes less fragmentation
			isactive[lp]= 2;
			if(b < high_ind) isactive[lp]= 1;
			asign(lp, b);
	//		cout << "LP index " << lp << " last-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}
		if(frag1 == frag2){				// Tie between first and last fit
			minIndAllocLP(lp);
	//		cout << "LP index " << lp << " minmax-asigned to index  " << spec_ind[lp] << endl;
			return 1;
		}
	}

//	tempAsign(lp, a);
	return 0;
}

int midFitAllocLP(int lp)
{
	int a, b, c;
	int p,q;

	int lp1;

	c = lp_size[lp];
	if(!c) return 0;

	int index1 = (high_ind + low_ind - c + 2) / 2;                 // Center index

	a= checkFirstFit(lp);				// Available index for first or last fit
	b= checkLastFit(lp);

	if(a == INF && b == INF){			// No available index to fit in
		blocked ++;
		return 0;
	}

	if(low_ind - high_ind < 4*req_Max+1){			//Avoiding blocking caused by middle fit LPs
		minFragAllocLP(lp);
	//	cout << "LP index " << lp << " minmax-asigned to index  " << spec_ind[lp] << endl;
		return 1;
	}else{										// Middle-fit allocation scheme
//		lpNode 	*cur = midFitList;
//		if( cur == NULL )									// Midlle is empty
		if( midlp == M){
			asign(lp, index1);								// Allocate to exact middle
			isactive[lp] = 3 ;
			//addToList(2, lp);
			midlp = lp;
//			cout << "LP index " << lp << " midfit-asigned0 to index  " << spec_ind[lp] << endl;
			return 1;
		}else{
//			cout << "LP index " << lp << ", a " << a << ",and frag1 , b " << b << ",and frag2 " << endl;
			double frag1 = 1;				//Firstfit
			double frag2 = 1;				//Lastfit

//			if(a < spec_ind[midlp])
				frag1 = checkFrag(lp, a);		// To avoid the rare cases where first(resp last) fit
//			if(b > spec_ind[midlp])
				frag2 = checkFrag(lp, b);		// need to allocate after the midlp

//			while( cur != NULL ){
//				lp1 = cur->x;
//				cur = cur->next;
				lp1 = midlp;
				double minfrag = frag1;
				if(frag2 < minfrag) minfrag = frag2;

				p = checkMoveUp(lp1);
			//	cout << "p: " << p << endl;
				q = checkMoveDown(lp1);
				double frag3 = checkFrag(lp1, p);						//Middle + move up
				if(frag3 < minfrag) minfrag = frag3;
				double frag4 = checkFrag(lp1, q);						//Middle + move down
				if(frag4 < minfrag) minfrag = frag4;

			//	cout << "Frags 1: " << frag1 << ", 2: " << frag2 << ", 3: " << frag3 << ", 4: " << frag4 << endl;

				if(minfrag == frag1){
					asign(lp, a);
					isactive[lp] = 1;
					if(a > index1) isactive[lp] = 2;
//					cout << "LP index " << lp << " firstfit-asigned to index  " << spec_ind[lp] << endl;
					return 1;
				}
				if(minfrag == frag2){
					asign(lp, b);
					isactive[lp] = 2;
					if(b < index1) isactive[lp] = 1;
//					cout << "LP index " << lp << " lastfit-asigned to index  " << spec_ind[lp] << endl;
					return 1;
				}
				if(minfrag == frag3){
					asign(lp, spec_ind[lp1]-c);
					isactive[lp] = 3;
					//addToList(2, lp);
					midlp = lp;
					isactive[lp1] = 2;
					//delFromList(2, lp1);
//					cout << "LP index " << lp << " midfit-asigned1 to index  " << spec_ind[lp] << endl;
					return 1;
				}
				if(minfrag == frag4){
					asign(lp, spec_ind[lp1]+lp_size[lp1]);
					isactive[lp] = 3;
				//	addToList(2, lp);
					midlp = lp;
					isactive[lp1] = 1;
				//	delFromList(2, lp1);
//					cout << "LP index " << lp << " midfit-asigned2 to index  " << spec_ind[lp] << endl;
					return 1;
				}

		//	}
		}
	}
	return 0;
}

int rpAllocLP(int lp)
{
	int a, b, c;

	int s = source[lp], d= dest[lp];

	c = lp_size[lp];
	if(!c) return 0;

	if(part[s][d] == 0){
		a= checkFirstFit(lp);				// Available index for first or last fit
		if(a == INF){			// No available index to fit in
			blocked ++;
			return 0;
		}
		if(a < low_ind){
	//	if(a != INF){
			isactive[lp]= 1;
			asign(lp, a);
//			cout << "LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}else{
			b= checkLastFit(lp);
			isactive[lp] = 2;
			asign(lp, b);
			return 1;
		}
//		if(lp == 561)
//			cout << "1: LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
	}else{
		b= checkLastFit(lp);
		if(b == INF){			// No available index to fit in
			blocked ++;
			return 0;
		}
		if(b > high_ind){
	//	if(b != INF){
			isactive[lp]= 2;
			asign(lp, b);
//			cout << "LP index " << lp << " last-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}else{
			a= checkFirstFit(lp);
			isactive[lp] = 1;
			asign(lp, a);
			return 1;
		}
	}
}

int firstFit1(int lp)
{
	int a=0, b=0, index=0;
	int i,j,p;
	bool asigned = 0, nofit = 0;
	int s,d;

	s = source[lp];
	d= dest[lp];
	b= lp_size[lp];
	index = 0;


	if(!b) return 0;

	a= checkExactFit(lp);
	if(a!=S){
		asign(lp, a);
		if(a <= high_ind+1) isactive[lp] = 1;
		if(a >= low_ind-1) isactive[lp] = 2;
		return 1;
	}

	while(index <= (S-b) && !asigned)   		  //Checking all spectrum range
	{
	//	cout << "LP checking index: " << lp << ", "<< index << endl;

		for(i=0;i<L;i++){							//Check path availibility
			if(spec[index][i] && path[s][d][i]){
				index++;
				break;
			}else{
				nofit = 0;
				for(j=0;j<b;j++){							//Checking if it fit (size)
					if(spec[index+j][i] && path[s][d][i]){
						index += j;
						nofit = 1;
						break;
					}
				}
			}
			if(nofit) break;
		}
		if(i==L && !nofit){
			isactive[lp]= 1;			// Side is defined before asign
			asign(lp, index);
			asigned= 1;
	//		cout << "LP index " << lp << " first-asigned to index" << spec_ind[lp] << endl;
		//	cout << "LP index is allocated #" << isactive[lp] << endl;
		}
	}
	if(!asigned){
		blocked ++;
		return 0;
	}else return 1;
}

int minIndAllocLP1(int lp)
{
	int a, b, c;

	c = lp_size[lp];
	if(!c) return 0;

	a= checkExactFit(lp);
	if(a!=S){
		asign(lp, a);
		if(a <= high_ind+1) isactive[lp] = 1;
		if(a >= low_ind-1) isactive[lp] = 2;
		return 1;
	}

	a= checkFirstFit(lp);				// Available index for first or last fit
	b= checkLastFit(lp);

	if(a == INF && b == INF){			// No available index to fit in
		blocked ++;
		return 0;
	}

	{	// min max index Allocation policy
		if(b > low_ind && a+c < high_ind){			// No added index in both side
			if(high_ind > S-low_ind){ 		// Lower side higher
				isactive[lp]= 2;
				asign(lp, b);
	//			cout << "LP index " << lp << " last-asigned_0 to index  " << spec_ind[lp] << endl;
				return 1;
			}else{							// Upper side higher
				isactive[lp]= 1;
				asign(lp, a);
	//			cout << "LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
				return 1;
			}
		}else{							// Added index in at least one side
			if(a + c -1 -high_ind < low_ind - b){ 		// Added index lower for first fit
				isactive[lp]= 1;
				asign(lp, a);
	//			cout << "LP index " << lp << " first-asigned_1 to index " << spec_ind[lp] << endl;
				return 1;
			}
			if(a + c -1 -high_ind > low_ind - b){ 		// Added index lower for last fit
				isactive[lp]= 2;
				asign(lp, b);
	//			cout << "LP index " << lp << " last-asigned_1 to index " << spec_ind[lp] << endl;
				return 1;
			}
			if(a + c -1 - high_ind == low_ind - b && high_ind +1 > S-low_ind){ 		// Added same for first and last fit, but last fit shorter
				isactive[lp]= 2;
				asign(lp, b);
	//			cout << "LP index " << lp << " last-asigned_2 to index " << spec_ind[lp] << endl;
				return 1;
			}else{										// first-fit shorter or the same size
				isactive[lp]= 1;
				asign(lp, a);
	//			cout << "LP index " << lp << " first-asigned_2 to index " << spec_ind[lp] << endl;
				return 1;
			}
	//	cout << "isactive a=" << isactive[lp] << endl;
		}
	}
}

int minFragAllocLP1(int lp)
{
	int a, b, c;

	c = lp_size[lp];
	if(!c) return 0;

	a= checkExactFit(lp);
	if(a!=S){
		asign(lp, a);
		if(a <= high_ind+1) isactive[lp] = 1;
		if(a >= low_ind-1) isactive[lp] = 2;
		return 1;
	}

	a= checkFirstFit(lp);				// Available index for first or last fit
	b= checkLastFit(lp);

	if(a == INF && b == INF){			// No available index to fit in
		blocked ++;
		return 0;
	}

	{									// min max index Allocation policy
		double frag1 = checkFrag(lp, a);
		double frag2 = checkFrag(lp, b);
	//	cout << "LP index " << lp << ", a " << a << ",and frag1 " << frag1<< ", b " << b << ",and frag2 " << frag2 << endl;

		if(frag1 < frag2){				// First fit causes less fragmentation
			isactive[lp]= 1;
			if(a >= low_ind) isactive[lp]= 2;
			asign(lp, a);
	//		cout << "LP index " << lp << " first-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}
		if(frag2 < frag1){				// Last fit causes less fragmentation
			isactive[lp]= 2;
			if(b < high_ind) isactive[lp]= 1;
			asign(lp, b);
	//		cout << "LP index " << lp << " last-asigned_0 to index  " << spec_ind[lp] << endl;
			return 1;
		}
		if(frag1 == frag2){				// Tie between first and last fit
			minIndAllocLP(lp);
	//		cout << "LP index " << lp << " minmax-asigned to index  " << spec_ind[lp] << endl;
			return 1;
		}
	}

//	tempAsign(lp, a);
	return 0;
}

int midFitAllocLP1(int lp)
{
	int a, b, c;
	int p,q;

	int lp1;

	c = lp_size[lp];
	if(!c) return 0;

	int index1 = (high_ind + low_ind - c + 2) / 2;                 // Center index

	if(low_ind - high_ind < 4*req_Max+1){			//Avoiding blocking caused by middle fit LPs
		minFragAllocLP(lp);
	//	cout << "LP index " << lp << " minmax-asigned to index  " << spec_ind[lp] << endl;
		return 1;
	}else{		// Middle-fit allocation scheme
	    a= checkExactFit(lp);
		if(a!=S){
			asign(lp, a);
			if(a <= high_ind+1) isactive[lp] = 1;
			if(a >= low_ind-1) isactive[lp] = 2;
			return 1;
		}

		a= checkFirstFit(lp);				// Available index for first or last fit
		b= checkLastFit(lp);

		if(a == INF && b == INF){			// No available index to fit in
			blocked ++;
			return 0;
		}

//		lpNode 	*cur = midFitList;
//		if( cur == NULL )									// Midlle is empty
		if( midlp == M){
			asign(lp, index1);								// Allocate to exact middle
			isactive[lp] = 3 ;
			//addToList(2, lp);
			midlp = lp;
//			cout << "LP index " << lp << " midfit-asigned0 to index  " << spec_ind[lp] << endl;
			return 1;
		}else{
//			cout << "LP index " << lp << ", a " << a << ",and frag1 , b " << b << ",and frag2 " << endl;
			double frag1 = 1;				//Firstfit
			double frag2 = 1;				//Lastfit

//			if(a < spec_ind[midlp])
				frag1 = checkFrag(lp, a);		// To avoid the rare cases where first(resp last) fit
//			if(b > spec_ind[midlp])
				frag2 = checkFrag(lp, b);		// need to allocate after the midlp

//			while( cur != NULL ){
//				lp1 = cur->x;
//				cur = cur->next;
				lp1 = midlp;
				double minfrag = frag1;
				if(frag2 < minfrag) minfrag = frag2;

				p = checkMoveUp(lp1);
			//	cout << "p: " << p << endl;
				q = checkMoveDown(lp1);
				double frag3 = checkFrag(lp1, p);						//Middle + move up
				if(frag3 < minfrag) minfrag = frag3;
				double frag4 = checkFrag(lp1, q);						//Middle + move down
				if(frag4 < minfrag) minfrag = frag4;

			//	cout << "Frags 1: " << frag1 << ", 2: " << frag2 << ", 3: " << frag3 << ", 4: " << frag4 << endl;

				if(minfrag == frag1){
					asign(lp, a);
					isactive[lp] = 1;
					if(a > index1) isactive[lp] = 2;
//					cout << "LP index " << lp << " firstfit-asigned to index  " << spec_ind[lp] << endl;
					return 1;
				}
				if(minfrag == frag2){
					asign(lp, b);
					isactive[lp] = 2;
					if(b < index1) isactive[lp] = 1;
//					cout << "LP index " << lp << " lastfit-asigned to index  " << spec_ind[lp] << endl;
					return 1;
				}
				if(minfrag == frag3){
					asign(lp, spec_ind[lp1]-c);
					isactive[lp] = 3;
					//addToList(2, lp);
					midlp = lp;
					isactive[lp1] = 2;
					//delFromList(2, lp1);
//					cout << "LP index " << lp << " midfit-asigned1 to index  " << spec_ind[lp] << endl;
					return 1;
				}
				if(minfrag == frag4){
					asign(lp, spec_ind[lp1]+lp_size[lp1]);
					isactive[lp] = 3;
				//	addToList(2, lp);
					midlp = lp;
					isactive[lp1] = 1;
				//	delFromList(2, lp1);
//					cout << "LP index " << lp << " midfit-asigned2 to index  " << spec_ind[lp] << endl;
					return 1;
				}

		//	}
		}
	}
	return 0;
}

int printSpec()
{
	int i,j;
	cout << "1:high_ind " << high_ind << ", low_ind " << low_ind << endl;
	cout << " Spectrum :=" << endl;						//Just for checking
	cout << "  l :";
	for (i=0;i<L;i++) cout <<"  "<< i ;
	cout << endl;
	for (i=S-1; i>=0; i--){
		if(i / 10 < 1) cout << " ";
		if(i / 100 <1) cout << " ";
		cout << i << " :";
		for(j=0;j<L;j++){
			cout << "  " << spec[i][j];
		}
		cout << endl;
	}
	cout << "  l :";
	for (i=0;i<L;i++) cout <<"  "<< i ;
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
		for(p=0;p<L;p++){
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
		for(p=0;p<L;p++){
			if(spec[index+j][p] == 1 && path_rr[p][lp] ==1) throw "プライマリパス割り当てエラー";
			spec[index+j][p] = spec[index+j][p] || path_rr[p][lp];
		}
	}

	return 0;
}

double checkFrag(int lp, int index)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;

	bool temp_spec[S][L];
	double pathFrag[N][N];
	double specFrag =0;
	double totalfrag=0, avefrag=0;
	int totalhops=0;

	for(i=0;i<N;i++){				// Initialise
		for(j=0;j<N;j++) 	pathFrag[i][j] = 0;
	}

	for(i=0;i<S;i++){				// Duplicating spectrum
		for(p=0;p<L;p++)
			temp_spec[i][p] = spec[i][p];
	}

	for(j=0;j<b;j++){				// Temporary asigning
		for(p=0;p<L;p++) temp_spec[index+j][p] = temp_spec[index+j][p] || path[s][d][p];
	}

	// cout << " Spectrum :=" << endl;						//Just for checking
	// cout << "l:";
	// for (i=0;i<=L;i++) cout <<"  "<< i ;
	// cout << endl;
	// for (i=S-1 ; i>=0; i--){
		// cout << i << " :";
		// for(j=0;j<L;j++){
			// cout << "  " << temp_spec[i][j];
		// }
		// cout << endl;
	// }
	// cout << endl;

	// cout << "1:high_ind " << high_ind << ", low_ind " << low_ind << endl;

	for(j=high_ind+b+2;j<low_ind-b;j++){				// Temporary removal of midfit LPs
		for(p=0;p<L;p++) temp_spec[j][p] = 0;
	}

	// cout << " Spectrum :=" << endl;						//Just for checking
	// cout << "l:";
	// for (i=0;i<=L;i++) cout <<"  "<< i ;
	// cout << endl;
	// for (i=S-1 ; i>=0; i--){
		// cout << i << " :";
		// for(j=0;j<L;j++){
			// cout << "  " << temp_spec[i][j];
		// }
		// cout << endl;
	// }
	// cout << endl;

	for(s1=0;s1<N;s1++){					// Determine avalaible SB and max SB
		for(d1=0;d1<N;d1++){						// for all paths
			int maxSB = 0, avSB = 0;
			int nonalign = 0, consecSB =0, maxcons=0;

//			if(linked_path[s][d][s1][d1]){			// For shared path
				for(i=0;i<S;i++){					// Checking spectrum for available aligned SB
					nonalign =0 ;
					for(j=0;j<L;j++){
						if(path[s1][d1][j]){		// Path s1-d1 using link j
							if(!temp_spec[i][j])  avSB++;		// Available SB i on link j
							if(temp_spec[i][j])  nonalign =1;	// SB i not aligned through the path of s1-d1
						}
					}
					if(nonalign==0) consecSB++;		// Increasing consecutive aligned SB
					if(nonalign || i==S-1){			// End of consecutive aligned SB
						if(maxcons < consecSB){
							maxcons = consecSB;
						}
						consecSB = 0;
					}
				}
				maxSB = hops[s1][d1]*maxcons;
				if(!avSB) pathFrag[s1][d1] = 0;
				if(avSB) pathFrag[s1][d1]= 1 - double(maxSB) / avSB;
				if(lp==8){
			//		cout << "Hops[" << s1 << "][" <<d1<<"]= " << hops[s1][d1] <<" maxSB and avSB: "<< maxSB <<" " <<avSB  << endl;
			//		cout << "1: pathFrag[" << s1 << "][" <<d1<<"]= " << pathFrag[s1][d1] << endl;
				}
//			}
		}
	}

	for(i=0;i<N;i++){				// Determining fragmentation
		for(j=0;j<N;j++){
			if(linked_path[s][d][i][j]){
//				cout << "2: pathFrag[" << i << "][" <<j<<"]= " << pathFrag[i][j] << endl;
				totalfrag += hops[i][j]*pathFrag[i][j];
				totalhops += hops[i][j];
			}
		}
	}
	specFrag = totalfrag / totalhops;
//		cout << "specFrag= " << specFrag << endl;

	return specFrag;
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

	for(i=0;i<S;i++){		// Checking spectrum for available aligned SB
		nonalign =0 ;
		for(j=0;j<L;j++){//全てのリンクについて
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
		if(i==S-1){//もし最高帯域スロットに達したならば
			if(consecSB == b)//もし必要な帯域スロットと同じならば
				return i-b+1;
			consecSB = 0;
		}
	}

	return S;
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

	for(i=0;i<S;i++){		// Checking spectrum for available aligned SB
		nonalign =0 ;
		for(j=0;j<L;j++){//全てのリンクについて
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
		if(i==S-1){//もし最高帯域スロットに達したならば
			if(consecSB == b)//もし必要な帯域スロットと同じならば
				return i-b+1;
			consecSB = 0;
		}
	}

	return S;
}

int checkExactFitRerouting(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;
	bool isGetRoot;

	for(i=0;i<S-b;i++){
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
	return S;
}

int checkExactBpRerouting(int lp)
{
	int i,j,p;
	int s = source[lp], d= dest[lp], b= lp_size[lp];
	int lp1, index1, s1, d1;
	bool isGetRoot;

	for(i=0;i<S-b;i++){
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
	return S;
}

int getPrimRoot(int s, int lp)
{
	int i, j, k;
	bool unconnectedFlag;
	int a, b;
	int to, cost;

	priority_queue<Node, vector<Node>, greater<Node> > q; // 優先度付き待ち行列

	Node Nodes[N];
	Node doneNode;

	//ノード情報input
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			if(link[i][j]>L) continue;
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
	for(i=0; i<N; i++){
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

	int path_rr_prev[N]; //to compare the routes

	if (b < N) {
		for (j = 0; j < N; j++) {//initialize
			for (k = 0; k < N; k++) {
				if(link[j][k]>L)continue;
				path_rr_prev[link[j][k]] = path_rr[link[j][k]][lp];
				path_rr[link[j][k]][lp] = 0;
			}
		}
		bool isSame = true;
		while (b < N && a != b) {
			path_rr[link[b][a]][lp] = 1;
			if (path_rr_prev[link[b][a]] != 1)
            {
                    isSame = false;
            }
			std::cout << "path_rr[" << link[b][a] << "][" << lp << "] = " << path_rr[link[b][a]][lp] << '\n';
			a = b;
			b = Nodes[a].from;
		}
		if (isSame)
		{
			rerouteOp++;
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

	Node Nodes[N];
	Node doneNode;

	//ノード情報input
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			if(link[i][j]>L) continue;
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

	for(i= 0; i < N ; i++){
	    for(j=0;j < Nodes[i].edges_to.size() ; j++){
	      // cout << "Nodes[" << i << "] = " << Nodes[i].edges_to[j] << '\n';
	    }
 	}

	// 初期化
	for(i=0; i<N; i++){
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

	if (b < N && a != b) {
		for (j = 0; j < N; j++) {//初期化
			for (k = 0; k < N; k++) {
				if(link[j][k]>L)continue;
				bp_rr[link[j][k]][lp] = 0;
			}
		}
		while (b < N && a != b) {
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

int removeLP_0(int lp)
{
	int s = source[lp], d= dest[lp];
	int index = spec_ind[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];

//	cout << "loop0, a=" << a << endl;

	if(a){										// a xor a =0, therefore if a lp is active
		for(i=0; i<b; i++){									// an xor with its path will remove it
			for(j=0;j<L;j++){
				spec[index+i][j] = path[s][d][j] ^	spec[index+i][j];
			}
		}

	//	cout << "loop2" << endl;

	isactive[lp] = 0;
	if(lp==midlp) midlp = M;
	//	cout << "S" << lp << " removed." << endl;
	//	printSpec();
	}
	return 0;
}

int removeLP(int lp)
{
	int s = source[lp], d= dest[lp];
	int index = spec_ind[lp],  b= lp_size[lp];
	int i,j;
	int a = isactive[lp];

	//		cout << "loop1, a=" << a << endl;

	if(a){										// a xor a =0, therefore if a lp is active
		for(i=0; i<b; i++){									// an xor with its path will remove it
			for(j=0;j<L;j++){
				spec[index+i][j] = path[s][d][j] ^	spec[index+i][j];
			}
		}

	//	cout << "loop2" << endl;

	isactive[lp] = 0;
	if(lp==midlp) midlp = M;
	//	cout << "S" << lp << " removed." << endl;
	//	printSpec();
	}

//	if(algoCall && algoCall != 2)
//		retune();
//	printSpec() ;
	return 0;
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
				for(j=0;j<L;j++){//全てのリンクに関して
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
				for(j=0;j<L;j++){
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
				for(j=0;j<L;j++){//全てのリンクに関して
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
				for(j=0;j<L;j++){
					if(spec[index+i][j] == 0 && bp_rr[j][lp] == 1) throw "バックアップパス消去エラー";
					spec[index+i][j] = bp_rr[j][lp] ^ spec[index+i][j];//pathが1ならspecは1
				}
			}
		}
	}

}

void fixMiddle()
{
	int s, d;
	int index, b;
	int i,j,p;
  	bool eol=0;

	if(midlp!=M && isactive[midlp]==3){
		index = spec_ind[midlp];
		b = lp_size[midlp];
		s = source[midlp];
		d = dest[midlp];

		for(i=0; i<b; i++){									// an xor with its path will remove it
			for(j=0;j<L;j++){
				spec[index+i][j] = path[s][d][j] ^	spec[index+i][j];
			}
		}

		int index1 = (high_ind + low_ind - b + 2) / 2;
		asign(midlp, index1);
	}

}

int checkMoveUp(int lp){
	int s, d;
	int index, a = lp_size[lp];
	int i,j,p;
  	bool eol=0;
	int k=0;

	i = lp;
	index = spec_ind[lp] + a;

	if(isactive[lp]==3){
		s = source[lp];
		d = dest[lp];
		eol = 0;

		while(index < S && !eol){		// Check to what extend it can be retuned
			for(j=0;j<L;j++){							//Check path availibility
				if(spec[index][j] && path[s][d][j]){
					eol = 1;
					break;
				}
			}
			if(eol) break;
			index++;
		}
	}

	return (index - a);

}

int checkMoveDown(int lp){
	int s, d;
	int index;
	int i,j,p;
  	bool eol=0;
	int k=0;

	index = spec_ind[lp]-1;

	i = lp;
	if(isactive[lp] == 3){
		s = source[lp];
		d = dest[lp];
		eol = 0;
		while(index >= 0 && !eol){
			for(j=0;j<L;j++){							//Check path availibility
				if(spec[index][j] && path[s][d][j]){
					eol = 1;
					break;
				}
			}
			if(eol) break;
			index--;
		}
	}

	return	(index + 1);
}
