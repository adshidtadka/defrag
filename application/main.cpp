#include "params.h"
#include "func.h"

using namespace std;

int main(int argc, char* argv[])
{
	int k;
	int a, b;
	int blockedff=0, blockedffh=0;

	int lp;

	// open output file
	ofstream ofs_result_txt;
	ofs_result_txt.open("./../result/result.txt");
	if(!ofs_result_txt){
		cout << "[error] Cannot open ./../result/result.txt" << endl;
		return 1;
	}

    ofstream ofs_result_csv;
    ofs_result_csv.open("./../result/result.csv");
    if(!ofs_result_csv){
            cout << "[error] Cannot open ./../result/result.csv" << endl;
            return 1;
    }

    ofstream operation_num_csv;
    operation_num_csv.open("./../result/operation_num.csv");
    if(!operation_num_csv){
            cout << "[error] Cannot open ./../result/operation_num.csv file" << endl;
            return 1;
    }

	ofs_result_txt 	<< "Simulation for NODE_NUM= "<< NODE_NUM <<", LINK_NUM = " << LINK_NUM << ", CAPASITY= "<< CAPASITY <<", REQ_SIZE_MAX=" << REQ_SIZE_MAX <<", DEFRAG_INTERVAL = "<< DEFRAG_INTERVAL <<", temp_max= "<< temp_max << ", DEFRAG_TIME = " << DEFRAG_TIME << ", MAX_HOP_NUM = " << MAX_HOP_NUM << endl;
	cout 			<< "Simulation for NODE_NUM= "<< NODE_NUM <<", LINK_NUM = " << LINK_NUM << ", CAPASITY= "<< CAPASITY <<", REQ_SIZE_MAX=" << REQ_SIZE_MAX <<", DEFRAG_INTERVAL = "<< DEFRAG_INTERVAL <<", temp_max= "<< temp_max << ", DEFRAG_TIME = " << DEFRAG_TIME << ", MAX_HOP_NUM = " << MAX_HOP_NUM << endl;
	ofs_result_txt 	<< endl << "REQ_NUM= "<< REQ_NUM << endl;
	cout 			<< endl << "REQ_NUM= "<< REQ_NUM << endl;


	for(int l=0; l<10; l++){
		int load = START_LOAD;
		if(l != 0){
			load = load + 10;
		}
		ofs_result_txt 	<< endl << "load= "<< load << endl;
		cout 			<< endl << "load= "<< load << endl;
		initialize();
		if (readInput(argc, &argv[0], load) == 1){
		    cout << "[error] cannot read input" << endl;
		    return 1;
		}
		int blkdItNone[10] = {}, blkdItConv[10] = {}, blkdItProp[10] = {}, blkdItConvIlp[10] = {}, blkdItPropIlp[10] = {};
		int togOpItConv[10] = {}, togOpItProp[10] = {};
		int realOpItConv[10] = {}, realOpItProp[10] = {};
		int rerouteOpItProp[10] = {};
		double stdItNoneArr[10][ITERATION] = {{}}, stdItConvArr[10][ITERATION] = {{}}, stdItPropArr[10][ITERATION] = {{}}, stdItConvIlpArr[10][ITERATION] = {{}}, stdItPropIlpArr[10][ITERATION] = {{}};
		for(int it = 0; it<ITERATION; it++){ //ランダムサンプルの数_5回				
				ofs_result_txt << endl << " Iteration= "<< it << endl;//何回めなのかを表示
				cout << endl << " Iteration= "<< it << endl;//何回めなのかを表示
				if(it==0) seed1 =  1448601515;  //time(NULL); 1444196111 1419542268 1424246601
				//最初は種を設定する
				seed1 = (it+1) * seed1;
				seed2 = (it+2) * seed1;
				temp_max = DEFRAG_TOTAL_TIME_MAX;//新しいパスが来ないときのデフラグ時間

				reInitialize();//経路とパス以外をゼロにする
				genDemands(load);//10万のパス情報を取り直す

				for(k=1; k<=1; k++){
					reInitialize();//経路とパスの情報をゼロにする

					initializeEvent(); //initialize startEvent endEvent defragEvent

						for(int i=0; i<REQ_NUM; i++){
							t_req[i]= double (t_req[i]);
							t_exp[i]= double (t_exp[i]);
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
						
						cout << "endEvent[" << REQ_NUM-1 << "].time = " << endEvent[REQ_NUM-1].time << ", DEFRAG_INTERVAL = " << DEFRAG_INTERVAL << endl;

						//make defrag event
						defragCount = round(endEvent[REQ_NUM-1].time/DEFRAG_INTERVAL);
						cout << "defragCount = " << defragCount << endl;
						defragEvent.clear();
						defragEvent.resize(defragCount);
						for (int c = 0; c < defragCount ; c++){
							defragEvent[c].time = c * DEFRAG_INTERVAL;
							defragEvent[c].type = 2;
						}
						// cout << "defragEvent[" << defragCount-1 << "].time = " << defragEvent[defragCount-1].time << endl;

					for(int j=0; j<=2; j++){	// To do all algorithms sequentially連続してアルゴリズムを実行する c=0なら2回ループ
						algoCall = j;

						//push to the priority queue
						// eventQueue.push(startEvent[0]);
						for (int i = 0; i < REQ_NUM; i++){
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
						// while(lp<REQ_NUM && t < T_end){//Start!!		// For running period 信号の数やステップ数が上限を超えない限り反復
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
											retuneBp(load);//デフラグを行う
										}
										catch(const char* err){
											cout << "ERR:ブロッキング中における" << err << endl;
											return 1;
										}
									}

									//新しくきたパスを評価する
									int s= source[nowEvent.lpNum], d= dest[nowEvent.lpNum], sort_val1 = 0, sort_val2 = 0;
									//ブロッキングが起きていなければパスをリストに加える
									if(b) addToList(nowEvent.lpNum, sort_val1, sort_val2);//activeList, backupList, mixtListに加えられる
									//最終パスならばシミュレーション終了
									if (nowEvent.lpNum == REQ_NUM-1){
										while(!eventQueue.empty()){
											eventQueue.pop();
										}
									}
	 							}
	 						}
							if(nowEvent.type == 2){//デフラグメンテーションがret_intごとに行われる
								try{
									retuneBp(load);//デフラグを行う
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
   						ofs_result_txt << "Simulation time = " << (double)(end - start) / CLOCKS_PER_SEC << " sec" <<endl;

						if(j==0){							// Reinitializing for next loop
							blockedff = blocked;
							blkdItNone[k] += blocked;
							stdItNoneArr[k][it] = blocked;
							cout << "Blocked request: First-fit w/o defragment                " << blockedff << endl;
							ofs_result_txt << "Blocked request: First-fit w/o defragment                " << blockedff << endl << endl;
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
							ofs_result_txt << "Blocked request: First-fit w/ stat defragment            " << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs_result_txt << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs_result_txt << "Number of reallocating operations:                       " << realOp << endl << endl;
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
							ofs_result_txt << "Blocked request: First-fit w/ stat rerouting defragment  " << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs_result_txt << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs_result_txt << "Number of reallocating operations:                       " << realOp << endl << endl;
							cout << "Number of rerouting operations:                          " << rerouteOp << endl;
                            ofs_result_txt << "Number of rerouting operations:                          " << rerouteOp << endl << endl;
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
							ofs_result_txt << "Blocked request: First-fit w/ stat ilp defragment        " << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs_result_txt << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs_result_txt << "Number of reallocating operations:                       " << realOp << endl << endl;
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
							ofs_result_txt << "Blocked request: First-fit w/ stat ilp reroute defragment" << blockedff << endl << endl;
							cout << "Number of switching operations:                          " << togOp << endl;
							ofs_result_txt << "Number of switching operations:                          " << togOp << endl;
							cout << "Number of reallocating operations:                       " << realOp << endl;
							ofs_result_txt << "Number of reallocating operations:                       " << realOp << endl << endl;
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
				ofs_result_csv << load << ",";
				operation_num_csv << load << ",";

				//信頼区間の計算
				double stdItNone, stdItConv, stdItProp, stdItConvIlp, stdItPropIlp;
				double stdItNoneDiff, stdItConvDiff, stdItPropDiff, stdItConvIlpDiff, stdItPropIlpDiff;
				stdItNone = standard(stdItNoneArr[k], ITERATION);
				stdItConv = standard(stdItConvArr[k], ITERATION);
				stdItProp = standard(stdItPropArr[k], ITERATION);
				stdItConvIlp = standard(stdItConvIlpArr[k], ITERATION);
				stdItPropIlp = standard(stdItPropIlpArr[k], ITERATION);
				stdItNoneDiff = blkdItNone[k]*0.05 - stdItNone*1.96;
				stdItConvDiff = blkdItConv[k]*0.05 - stdItConv*1.96;
				stdItPropDiff = blkdItProp[k]*0.05 - stdItProp*1.96;
				stdItConvIlpDiff = blkdItConvIlp[k]*0.05 - stdItConvIlp*1.96;
				stdItPropIlpDiff = blkdItPropIlp[k]*0.05 - stdItPropIlp*1.96;

				cout << "Av blocked request blkdItNone:         " << blkdItNone[k]/ITERATION << endl;
				ofs_result_txt << "Av blocked request blkdItNone:         " << blkdItNone[k]/ITERATION << endl;
				ofs_result_csv << blkdItNone[k]/ITERATION << ",";
				cout << "Confidence interval blkdItNone:        " << stdItNone << endl;
				ofs_result_txt << "Confidence interval blkdItNone:        " << stdItNone << endl;
				cout << "Av blocked request blkdItConv:         " << blkdItConv[k]/ITERATION << endl;
				ofs_result_txt << "Av blocked request blkdItConv:         " << blkdItConv[k]/ITERATION << endl;
				ofs_result_csv << blkdItConv[k]/ITERATION << ",";
				cout << "Confidence interval blkdItConv:        " << stdItConv << endl;
				ofs_result_txt << "Confidence interval blkdItConv:        " << stdItConv << endl;
				cout << "Av blocked request blkdItProp:         " << blkdItProp[k]/ITERATION << endl;
				ofs_result_txt << "Av blocked request blkdItProp:         " << blkdItProp[k]/ITERATION << endl;
				ofs_result_csv << blkdItProp[k]/ITERATION << ",";
				cout << "Confidence interval blkdItProp:        " << stdItProp << endl;
				ofs_result_txt << "Confidence interval blkdItProp:        " << stdItProp << endl;
				cout << "Av blocked request blkdItConvIlp:      " << blkdItConvIlp[k]/ITERATION << endl;
				ofs_result_txt << "Av blocked request blkdItConvIlp:      " << blkdItConvIlp[k]/ITERATION << endl;
				ofs_result_csv << blkdItConvIlp[k]/ITERATION << ",";
				cout << "Confidence interval blkdItConvIlp:     " << stdItConvIlp << endl;
				ofs_result_txt << "Confidence interval blkdItConvIlp:     " << stdItConvIlp << endl;
				cout << "Av blocked request blkdItPropIlp:      " << blkdItPropIlp[k]/ITERATION << endl;
				ofs_result_txt << "Av blocked request blkdItPropIlp:      " << blkdItPropIlp[k]/ITERATION << endl;
				ofs_result_csv << blkdItPropIlp[k]/ITERATION << endl;
				cout << "Confidence interval blkdItPropIlp:     " << stdItPropIlp << endl;
				ofs_result_txt << "Confidence interval blkdItPropIlp:     " << stdItPropIlp << endl;
				cout << "Av toggle operations togOpItConv:      " << togOpItConv[k]/ITERATION << endl;
				ofs_result_txt << "Av toggle operations togOpItConv:      " << togOpItConv[k]/ITERATION << endl;
				operation_num_csv << togOpItConv[k]/ITERATION << ",";
				cout << "Av move operations realOpItConv:       " << realOpItConv[k]/ITERATION << endl;
				ofs_result_txt << "Av move operations realOpItConv:       " << realOpItConv[k]/ITERATION << endl;
				operation_num_csv << realOpItConv[k]/ITERATION << ",";
				cout << "Av toggle operations togOpItProp:      " << togOpItProp[k]/ITERATION << endl;
				ofs_result_txt << "Av toggle operations togOpItProp:      " << togOpItProp[k]/ITERATION << endl;
				operation_num_csv << togOpItProp[k]/ITERATION << ",";
				cout << "Av move operations realOpItProp:       " << realOpItProp[k]/ITERATION << endl;
				ofs_result_txt << "Av move operations realOpItProp:       " << realOpItProp[k]/ITERATION << endl;
            	operation_num_csv << realOpItProp[k]/ITERATION << ",";
            	cout << "Av reroute operations rerouteOpItProp: " << rerouteOpItProp[k]/ITERATION << endl << endl;
            	ofs_result_txt << "Av reroute operations rerouteOpItProp: " << rerouteOpItProp[k]/ITERATION << endl << endl;
            	operation_num_csv << rerouteOpItProp[k]/ITERATION << endl;
            }
		}
	ofs_result_txt.close();

	lpNode *cur;
	cur = realList;
	while (cur != NULL) {
		deleteLP(cur->x, 0);
		cur = cur->next;
	}

	cout << endl;
	cout << "t: " << t << " Seed1: " << seed1 << ", Seed2: " << seed2 << endl<< endl;

	return 0;
}

int readInput(int argc, char* argv[0], int load)
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
