/* Add retuning time for 1+1 */
/* Stop retuning if new request */
/* Use speed up retuning */
/* Use T_temp, put high value if not needed */
/* include sorted list*/
/* Including backup first */
/* Parallel realloc used for conventional, Check retuning function*/

#define INF 999999999
#define T_end 1245201000        //28270000000

#define IT 3             // Number of ramdom sample

#define M 10000	     // Maximum number of demands

#define N 11          	// N number of Nodes.
#define L 52			// L number of directed Links. 有方向グラフ (N,L) = (11, 28), (5, 12), (14,44), (11, 52), (14, 46), (25, 84)
#define S 400			// S number of spec slots per link

#define A1 350		// Traffic load
#define H  10			// 1/mu, Average holding time in Tu
#define req_Max 16		// Maximum demand size 占有帯域スロット数

#define K 1000			//to explore less possibility of repeating of next arrival time of requests.  

#define DEFRAG_TIME 0.1			//to explore less possibility of repeating of next arrival time of requests.  

#define  T_temp 100     	// Vaiting time allowed to retune before adding new request
#define R_int 0.2		// Retuning period
// #define Sorting 1		// 0 if not used, 1 by length, 2 by size, 3 by blocks
//ソート方法が3種類与えられている
#define maxStep 3

#define spfact 1

#define LIMIT_HOP_NUM 99999

#include "allfunction_206.h"

using namespace std; //名前空間の指定

	// Declaring functions

int main(int argc, char* argv[])
{
	int i, j, k, l, p;
	int a, b;
	int blockedff=0, blockedffh=0;

	int lp;

	double ret_int = R_int;
	int Sorting;

	int it;

	lpNode *cur; //curというポインタ宣言, lpNodeは構造体
	//xはスロット番号, yは占有帯域の大きさ(ソート用),zはプライマリ・バックアップ判別子

	ofstream ofs1; //出力用ストリーム
	ofs1.open("./../result/blocked_205_b.txt"); //出力ファイルをオープンする
	if(!ofs1){
		cout<< "Cannot open blocked file"<<endl;
		return 1;
	}

    ofstream ofs4;
    ofs4.open("./../result/blocked_205_b.csv"); //blocking probability
    if(!ofs4){
            cout<< "Cannot open csv blocked file"<<endl;
            return 1;
    }

    ofstream ofs5;
    ofs5.open("./../result/operation_num.csv"); //operation number
    if(!ofs5){
            cout<< "Cannot open csv operation_num file"<<endl;
            return 1;
    }

	ofs1 << "Simulation for N= "<< N <<", Links = " << L << ", Spec= "<< S <<", req_Max=" << req_Max <<", A1= "<< A<< ", R_int = "<< R_int <<", temp_max= "<< temp_max << ", DEFRAG_TIME = " << DEFRAG_TIME << ", LIMIT_HOP_NUM = " << LIMIT_HOP_NUM << endl;
	cout << "Simulation for N= "<< N <<", Links = " << L << ", Spec= "<< S <<", req_Max=" << req_Max <<", A1= "<< A<< ", R_int = "<< R_int <<", temp_max= "<< temp_max << ", DEFRAG_TIME = " << DEFRAG_TIME << ", LIMIT_HOP_NUM = " << LIMIT_HOP_NUM << endl;
		//number of nodes, number of slots, maxi demand size, Traffic load, retuning period, Vaiting time allowed to retune before adding new request
	ofs1 << endl << " M= "<< M << endl;
		//max number of demand

	for(p=1; p<2; p++){ //Listを作る方法の違い
		Sorting= p;
		ofs1 << endl << "Sorting= "<< Sorting << endl;
		cout << endl;
		cout << endl << "Sorting= "<< Sorting << endl;

		A = A1;//現在の通信量
		ofs1 << endl << " A= "<< A << endl;
		cout << endl << " A= "<< A << endl;
		for(l=0; l<10; l++){ //通信量を変更するためのループ
			cout << "TERM NUMBER = " << l << endl;
			if(l){
				A = A + 10;//通信量
				ofs1 << endl << " A= "<< A << endl;
				cout << endl << " A= "<< A << endl;
			} //l=0のときは既にAを表示しているため実行しない
			initialize(); //色々全部0にする
			if (readInput(argc, &argv[0]) == 1){
                cout << "cannot read input" << endl;
                return 1;
            }//入力する情報を持って来る		
			int blkdItNone[10] = {}, blkdItConv[10] = {}, blkdItProp[10] = {}, blkdItConvIlp[10] = {}, blkdItPropIlp[10] = {};
			int togOpItConv[10] = {}, togOpItProp[10] = {};
			int realOpItConv[10] = {}, realOpItProp[10] = {};
			int rerouteOpItProp[10] = {};
			double stdItNoneArr[10][IT] = {{}}, stdItConvArr[10][IT] = {{}}, stdItPropArr[10][IT] = {{}}, stdItConvIlpArr[10][IT] = {{}}, stdItPropIlpArr[10][IT] = {{}};
			for(it = 0; it<IT; it++){ //ランダムサンプルの数_5回				
				ofs1 << endl << " Iteration= "<< it << endl;//何回めなのかを表示
				cout << endl << " Iteration= "<< it << endl;//何回めなのかを表示
				if(it==0) seed1 =  1448601515;  //time(NULL); 1444196111 1419542268 1424246601
				//最初は種を設定する
				seed1 = (it+1) * seed1;//種を変える
				seed2 = (it+2) * seed1;						//(k+2) * time(NULL);

				ret_int = R_int;//最大ステップ数
				temp_max = T_temp;//新しいパスが来ないときのデフラグ時間

				reInitialize();//経路とパス以外をゼロにする
				genDemands();//10万のパス情報を取り直す

				for(k=1; k<=1; k++){//k=1のみのループ
					reInitialize();//経路とパスの情報をゼロにする

					initializeEvent(); //initialize startEvent endEvent defragEvent

					// double spfact;//正規化のための変数
					// if(k==1) spfact = 0.01;
					// if(k==2) spfact = 2;
					// if(k==3) spfact = 2.5;
					// if(k==4) spfact = 2;

					if(k){
						ofs1 << endl << "Speeding = " << spfact << endl;
						cout << endl << "Speeding = " << spfact << endl;
						// startEvent[0].time = 0;
						// startEvent[0].type = 0;
						// startEvent[0].lpNum = 0;
						for(i=0; i<M; i++){//正規化
							t_req[i]= double (t_req[i]) * spfact;
							t_req_event[i+1] *= spfact;
							t_exp[i]= double (t_exp[i]) * spfact;
							t_exp_event[i] *= spfact;
							if(t_req[i] == t_exp[i]){
								t_exp[i]++; //到着時刻と切断時刻が同じであれば1ずらす
								// t_exp_event[i] += DBL_MIN;	
							}
							endEvent[i].time = t_exp_event[i];	//event driven
							endEvent[i].type = 1;
							endEvent[i].lpNum = i;
							// cout << "t_req_event[" << i << "] = " << t_req_event[i] << " t_exp_event[" << i << "] = " << t_exp_event[i] << endl;
							// cout << "startEvent[" << i+1 << "].time = " << startEvent[i+1].time << endl;
							startEvent[i].time = t_req_event[i];
							startEvent[i].type = 0;
							startEvent[i].lpNum = i;
						}

						ret_int = spfact * ret_int;
						temp_max = temp_max * spfact;
						
						cout << "endEvent[" << M-1 << "].time = " << endEvent[M-1].time << ", ret_int = " << ret_int << endl;
						//make defrag event
						defragCount = round(endEvent[M-1].time/ret_int);
						cout << "defragCount = " << defragCount << endl;
						defragEvent.clear();
						defragEvent.resize(defragCount);
						for (int c = 0; c < defragCount ; c++){
							defragEvent[c].time = c * ret_int;
							defragEvent[c].type = 2;
						}
						// cout << "defragEvent[" << defragCount-1 << "].time = " << defragEvent[defragCount-1].time << endl;
					}

					for(j=0; j<=2; j++){	// To do all algorithms sequentially連続してアルゴリズムを実行する c=0なら2回ループ
						algoCall = j;

						//push to the priority queue
						// eventQueue.push(startEvent[0]);
						for (i = 0; i < M; i++){
							eventQueue.push(endEvent[i]);
							// cout << "endEvent[" << i << "].time= " << endEvent[i].time << endl;
							eventQueue.push(startEvent[i]);
							// cout << "startEvent[" << i << "].time= " << startEvent[i].time << endl;
						}
						for (int c = 0; c < defragCount; c++){
							// cout << "defragEvent[" << c << "].time= " << defragEvent[c].time << endl;
							eventQueue.push(defragEvent[c]);
						}

						t = 0; lp=0;                  		// for time

						clock_t start, end; //time stump
						start = clock();

						/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						// while(lp<M && t < T_end){//Start!!		// For running period 信号の数やステップ数が上限を超えない限り反復
						while(!eventQueue.empty()){//Start!!		// For running period 信号の数やステップ数が上限を超えない限り反復
							b = 0;

							// if expired path is
							while(!deleteQueue.empty()) {
							    nowEvent = deleteQueue.top();
							    deleteQueue.pop();
							    try {
							    	removeLP1_1(nowEvent.lpNum, algoCall);
							    }
							    catch(const char* err) {
							    	cout << "ERR:パス切断中における" << err << endl;
							    	return 1;
							    }
								delFromList(1, nowEvent.lpNum);			
								delFromList(2, nowEvent.lpNum);			
							}

							// next event
							nowEvent = eventQueue.top();
	 						eventQueue.pop();
	 						t = nowEvent.time;
	 						// cout << "nowEvent.time = " << t << ", nowEvent.type = " << nowEvent.type << endl;

	 						if(nowEvent.type == 1){
	 							t = nowEvent.time;
	 							try {
							    	removeLP1_1(nowEvent.lpNum, algoCall);
							    }
							    catch(const char* err) {
							    	cout << "ERR:パス切断中における" << err << endl;
							    	return 1;
							    }
								delFromList(1, nowEvent.lpNum);							// Removing from linked list(active list)
								delFromList(2, nowEvent.lpNum);							// Removing from linked list(mixtlist)
	 						}
	 						if(nowEvent.type == 0){
	 							if(lp_size[nowEvent.lpNum]){
	 								last_lp = nowEvent.lpNum;//パスのIDをlast_lpに移す

	 								try{
										b = firstFit1_1(nowEvent.lpNum, algoCall);
									}
									catch(const char* err){
										cout << "ERR:パス割り当て中における" << err << endl;
										return 1;
									}

									if(!b){//もしブロッキングが起こってしまったら
										try{
											retuneBp();//デフラグを行う
										}
										catch(const char* err){
											cout << "ERR:ブロッキング中における" << err << endl;
											return 1;
										}
									}

									//新しくきたパスを評価する
									int s= source[nowEvent.lpNum], d= dest[nowEvent.lpNum], sort_val1 = 0, sort_val2 = 0;
									if(Sorting == 1){
										sort_val1 = hops[s][d];
										sort_val2 = bhops[s][d];
									}
									if(Sorting == 3){
										sort_val1 = lp_size[nowEvent.lpNum];
										sort_val2 = lp_size[nowEvent.lpNum];
									}
									if(Sorting == 2){
										sort_val1 = lp_size[nowEvent.lpNum] *hops[s][d];
										sort_val2 = lp_size[nowEvent.lpNum] *bhops[s][d];
									}
									//ブロッキングが起きていなければパスをリストに加える
									if(b) addToList(nowEvent.lpNum, sort_val1, sort_val2);//activeList, backupList, mixtListに加えられる
									//最終パスならばシミュレーション終了
									if (nowEvent.lpNum == M-1){
										while(!eventQueue.empty()){
											eventQueue.pop();
										}
									}
	 							}
	 						}
							if(nowEvent.type == 2){//デフラグメンテーションがret_intごとに行われる
								try{
									retuneBp();//デフラグを行う
								}
								catch(const char* err){
									cout << "ERR:デフラグ中における" << err << endl;
									return 1;
								}							
							}
						}//end!!!
						////////////////////////////////////////////////////////////////////////
						end = clock();
   						cout << "Simulation time = " << (double)(end - start) / CLOCKS_PER_SEC << " sec" <<endl;
   						ofs1 << "Simulation time = " << (double)(end - start) / CLOCKS_PER_SEC << " sec" <<endl;

						if(j==0){							// Reinitializing for next loop
							blockedff = blocked;
							blkdItNone[k] += blocked;
							stdItNoneArr[k][it] = blocked;
							cout << "Blocked request: First-fit w/o defragment                " << blockedff << endl;
							ofs1 << "Blocked request: First-fit w/o defragment                " << blockedff << endl << endl;
							// printSpec();
							// writeOutput();
							// statDefrag();
							// printSpec();
							reInitialize();
						}
						if(j==1){
							blockedff = blocked;
							blkdItConv[k] += blocked;
							stdItConvArr[k][it] = blocked;
							togOpItConv[k]  += togOp;
							realOpItConv[k] += realOp;
							cout << "Blocked request: First-fit w/ stat defragment            " << blockedff << endl;
							ofs1 << "Blocked request: First-fit w/ stat defragment            " << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs1 << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs1 << "Number of reallocating operations:                       " << realOp << endl << endl;
							// printSpec();
							// statDefrag();
							// printSpec();
							reInitialize();
						}
						if(j==2){
							blockedff = blocked;
							blkdItProp[k] += blocked;
							stdItPropArr[k][it] = blocked;
							togOpItProp[k]  += togOp;
							realOpItProp[k] += realOp;
							rerouteOpItProp[k] += rerouteOp;
							cout << "Blocked request: First-fit w/ stat rerouting defragment  " << blockedff << endl;
							ofs1 << "Blocked request: First-fit w/ stat rerouting defragment  " << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs1 << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs1 << "Number of reallocating operations:                       " << realOp << endl << endl;
							cout << "Number of rerouting operations:                          " << rerouteOp << endl;
                            ofs1 << "Number of rerouting operations:                          " << rerouteOp << endl << endl;
							// printSpec();
							// statDefrag();
							// printSpec();
							reInitialize();
						}
						if(j==3){
							blockedff = blocked;
							blkdItConvIlp[k] += blocked;
							stdItConvIlpArr[k][it] = blocked;
							cout << "Blocked request: First-fit w/ stat ilp defragment        " << blockedff << endl;
							ofs1 << "Blocked request: First-fit w/ stat ilp defragment        " << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs1 << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs1 << "Number of reallocating operations:                       " << realOp << endl << endl;
							// printSpec();
							// statDefrag();
							// printSpec();
							reInitialize();
						}
						if(j==4){
							blockedff = blocked;
							blkdItPropIlp[k] += blocked;
							stdItPropIlpArr[k][it] = blocked;
							cout << "Blocked request: First-fit w/ stat ilp reroute defragment" << blockedff << endl;
							ofs1 << "Blocked request: First-fit w/ stat ilp reroute defragment" << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs1 << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs1 << "Number of reallocating operations:                       " << realOp << endl << endl;
							// printSpec();
							// statDefrag();
							// printSpec();
							reInitialize();
						}
					}

				}

			}
			for(k=1; k<=1; k++){//k=1のみのループ

				//traffic load
				ofs4 << A << ",";
				ofs5 << A << ",";

				// double spfact;

				// if(k==1) spfact = 0.01;
				// if(k==2) spfact = 0.02;
				// if(k==3) spfact = 0.05;
				// if(k==4) spfact = 0.1;

				//speeding
				ofs4 << spfact << ",";
				ofs5 << spfact << ",";

				//信頼区間の計算
				double stdItNone, stdItConv, stdItProp, stdItConvIlp, stdItPropIlp;
				double stdItNoneDiff, stdItConvDiff, stdItPropDiff, stdItConvIlpDiff, stdItPropIlpDiff;
				stdItNone = standard(stdItNoneArr[k], IT);
				stdItConv = standard(stdItConvArr[k], IT);
				stdItProp = standard(stdItPropArr[k], IT);
				stdItConvIlp = standard(stdItConvIlpArr[k], IT);
				stdItPropIlp = standard(stdItPropIlpArr[k], IT);
				stdItNoneDiff = blkdItNone[k]*0.05 - stdItNone*1.96;
				stdItConvDiff = blkdItConv[k]*0.05 - stdItConv*1.96;
				stdItPropDiff = blkdItProp[k]*0.05 - stdItProp*1.96;
				stdItConvIlpDiff = blkdItConvIlp[k]*0.05 - stdItConvIlp*1.96;
				stdItPropIlpDiff = blkdItPropIlp[k]*0.05 - stdItPropIlp*1.96;

				cout << "Av blocked request blkdItNone:         " << blkdItNone[k]/IT << endl;
				ofs1 << "Av blocked request blkdItNone:         " << blkdItNone[k]/IT << endl;
				ofs4 << blkdItNone[k]/IT << ",";
				cout << "Confidence interval blkdItNone:        " << stdItNone << endl;
				ofs1 << "Confidence interval blkdItNone:        " << stdItNone << endl;
				cout << "Av blocked request blkdItConv:         " << blkdItConv[k]/IT << endl;
				ofs1 << "Av blocked request blkdItConv:         " << blkdItConv[k]/IT << endl;
				ofs4 << blkdItConv[k]/IT << ",";
				cout << "Confidence interval blkdItConv:        " << stdItConv << endl;
				ofs1 << "Confidence interval blkdItConv:        " << stdItConv << endl;
				cout << "Av blocked request blkdItProp:         " << blkdItProp[k]/IT << endl;
				ofs1 << "Av blocked request blkdItProp:         " << blkdItProp[k]/IT << endl;
				ofs4 << blkdItProp[k]/IT << ",";
				cout << "Confidence interval blkdItProp:        " << stdItProp << endl;
				ofs1 << "Confidence interval blkdItProp:        " << stdItProp << endl;
				cout << "Av blocked request blkdItConvIlp:      " << blkdItConvIlp[k]/IT << endl;
				ofs1 << "Av blocked request blkdItConvIlp:      " << blkdItConvIlp[k]/IT << endl;
				ofs4 << blkdItConvIlp[k]/IT << ",";
				cout << "Confidence interval blkdItConvIlp:     " << stdItConvIlp << endl;
				ofs1 << "Confidence interval blkdItConvIlp:     " << stdItConvIlp << endl;
				cout << "Av blocked request blkdItPropIlp:      " << blkdItPropIlp[k]/IT << endl;
				ofs1 << "Av blocked request blkdItPropIlp:      " << blkdItPropIlp[k]/IT << endl;
				ofs4 << blkdItPropIlp[k]/IT << endl;
				cout << "Confidence interval blkdItPropIlp:     " << stdItPropIlp << endl;
				ofs1 << "Confidence interval blkdItPropIlp:     " << stdItPropIlp << endl;
				cout << "Av toggle operations togOpItConv:      " << togOpItConv[k]/IT << endl;
				ofs1 << "Av toggle operations togOpItConv:      " << togOpItConv[k]/IT << endl;
				ofs5 << togOpItConv[k]/IT << ",";
				cout << "Av move operations realOpItConv:       " << realOpItConv[k]/IT << endl;
				ofs1 << "Av move operations realOpItConv:       " << realOpItConv[k]/IT << endl;
				ofs5 << realOpItConv[k]/IT << ",";
				cout << "Av toggle operations togOpItProp:      " << togOpItProp[k]/IT << endl;
				ofs1 << "Av toggle operations togOpItProp:      " << togOpItProp[k]/IT << endl;
				ofs5 << togOpItProp[k]/IT << ",";
				cout << "Av move operations realOpItProp:       " << realOpItProp[k]/IT << endl;
				ofs1 << "Av move operations realOpItProp:       " << realOpItProp[k]/IT << endl;
            	ofs5 << realOpItProp[k]/IT << ",";
            	cout << "Av reroute operations rerouteOpItProp: " << rerouteOpItProp[k]/IT << endl << endl;
            	ofs1 << "Av reroute operations rerouteOpItProp: " << rerouteOpItProp[k]/IT << endl << endl;
            	ofs5 << rerouteOpItProp[k]/IT << endl;
            }
		}
	}
	ofs1.close();

//	printSpec();

	cur = realList;
	while ( cur != NULL ) {
	//	if(bp_ind[cur->x] >= 16){
	//		cout << cur->x << " indexes " << spec_ind[cur->x] << ",  " << bp_ind[cur->x]<< endl;
	//		 printSpec();
	//	}
		deleteLP(cur->x, 0);
//		if(cur->z) cout<< cur->x << " prime index " << spec_ind[cur->x]<< endl;
//		if(!cur->z) cout<< cur->x << " backup index " << bp_ind[cur->x]<< endl;
		cur = cur->next;
	}

	cout << endl;

	cout << "t: " << t << " Seed1: " << seed1 << ", Seed2: " << seed2 << endl<< endl;

//	printSpec();
}

int readInput(int argc, char* argv[0])
{
	int i,j,k,l, p;
	int a, b;
	char tmp= 'r';

	ifstream fin; 									//reading from files　入力用ストリーム

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
	for(k=0; k<L; k++){
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

	for (i=0;i<N;i++){				//Comparing path, may be usefull
		for (j=0;j<N;j++){
			for (k=0;k<N;k++){
				for (l=0;l<N;l++){
					if(k == i && l == j) linked_path[i][j][k][l] = 0; //2組の発着ノードの組が同じ発着ノードであれば無視する
					else{
						for(p=0;p<L;p++){ //リンクの数だけforループを回す
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
	genDemands();
	return 0;
}

int writeOutput()
{
	int i,j,ind=0;
	int lp,n,f0;
	int s,d;
	int c = maxStep;
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
	ofs2 << "param M := " << m <<";" << endl;
	ofs2 << "param B := " << S <<";" << endl;
	ofs2 << "param L := " << L <<";" << endl;
	ofs2 << "param C := " << c <<";" << endl;
	ofs2 << endl;
	ofs2 << "param : S K N F0 := " << endl;

	ofstream ofs3;
	ofs3.open("running.txt");
	if(!ofs3){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}
	ofs3 << "Load A := " << A <<";" << endl;
	ofs3 << "param M := " << m <<";" << endl;
	ofs3 << "Last LP := " << lp <<";" << endl;
	ofs3.close();
	// cout << "Load A := " << A <<";" << endl;
	// cout << "param M := " << m <<";" << endl;
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

		for(j=0;j<L;j++){
			if(path[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
		for(j=0;j<L;j++){
			if(bp[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
	}
	ofs2 << ";" << endl << endl;

	ofs2.close();
	return 0;
}

int writeOutputPy()
{
	int i,j,ind=0;
	int lp,n,f0,k;
	int s,d;
	int c = maxStep;
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

	ofs2 << "param S   := " << m   << endl;
	ofs2 << "param absP:= " << m*2 << endl;
	ofs2 << "param absF:= " << S   << endl;
	ofs2 << "param absE:= " << L   << endl;
	ofs2 << "param T   := " << c   << endl;
	ofs2 << endl;

	ofstream ofs3;
	ofs3.open("./../result/ssr_lno_running.txt");
	if(!ofs3){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}
	ofs3 << "Load A := " << A <<";" << endl;
	ofs3 << "param M := " << m <<";" << endl;
	ofs3 << "Last LP := " << lp <<";" << endl;
	ofs3.close();
	// cout << "Load A := " << A <<";" << endl;
	// cout << "param M := " << m <<";" << endl;
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

		for(j=0;j<L;j++){
			if(path[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
		for(j=0;j<L;j++){
			if(bp[s][d][j]) ofs2 << ind <<" "<< j << endl;
		}
		ind++;
	}
	ofs2 << ";" << endl << endl;

	// printSpec();

	ofs2.close();
	return 0;
}

int writeOutputReroutingPy()
{
	int i,j,ind=0;
	int lp,n,f0,k;
	int s,d;
	int c = maxStep;
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

	ofs2 << "param S   := " << m   << endl;
	ofs2 << "param absP:= " << m*2 << endl;
	ofs2 << "param absF:= " << S   << endl;
	ofs2 << "param absE:= " << L   << endl;
	ofs2 << "param absV:= " << N   << endl;
	ofs2 << "param T   := " << c   << endl;
	ofs2 << endl;

	ofstream ofs3;
	ofs3.open("./../result/ssr_lno_running.txt");
	if(!ofs3){
		cout<< "Cannot open ouput file"<<endl;
		return 1;
	}
	ofs3 << "Load A := " << A <<";" << endl;
	ofs3 << "param M := " << m <<";" << endl;
	ofs3 << "Last LP := " << lp <<";" << endl;
	ofs3.close();
	// cout << "Load A := " << A <<";" << endl;
	// cout << "param M := " << m <<";" << endl;
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
	for ( i = 0; i < N; i++){
		for ( j = 0; j < N; j++){
			if (link[i][j] < L){
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
		for ( i = 0; i < N; i++) {
			for ( j = 0; j < N; j++ ){
				if(link[i][j] < L){
					if (path_rr[link[i][j]][lp]){
						ofs2 << ind << " " << i << " " << j << endl;
						// cout << ind << " " << i << " " << j << endl;
					}
				}
			}
		}

		ind++;//バックアップパスが連番になっている
		for ( i = 0; i < N; i++) {
			for ( j = 0; j < N; j++ ){
				if(link[i][j] < L){
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
