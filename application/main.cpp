#include "params.h"
#include "func.h"

int main(int argc, char* argv[])
{
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
            cout << "[error] Cannot open ./../result/operation_num.csv" << endl;
            return 1;
    }

	ofs_result_txt 	<< "Simulation for NODE_NUM= "<< NODE_NUM <<", LINK_NUM = " << LINK_NUM << ", CAPASITY= "<< CAPASITY << ", REQ_NUM= "<< REQ_NUM <<", REQ_SIZE_MAX=" << REQ_SIZE_MAX <<", DEFRAG_INTERVAL = "<< DEFRAG_INTERVAL <<", DEFRAG_TOTAL_TIME_MAX = "<< DEFRAG_TOTAL_TIME_MAX << ", DEFRAG_TIME = " << DEFRAG_TIME << ", MAX_HOP_NUM = " << MAX_HOP_NUM << endl;
	cout 			<< "Simulation for NODE_NUM= "<< NODE_NUM <<", LINK_NUM = " << LINK_NUM << ", CAPASITY= "<< CAPASITY << ", REQ_NUM= "<< REQ_NUM <<", REQ_SIZE_MAX=" << REQ_SIZE_MAX <<", DEFRAG_INTERVAL = "<< DEFRAG_INTERVAL <<", DEFRAG_TOTAL_TIME_MAX = "<< DEFRAG_TOTAL_TIME_MAX << ", DEFRAG_TIME = " << DEFRAG_TIME << ", MAX_HOP_NUM = " << MAX_HOP_NUM << endl;

	int load = LOAD_START;
	for(int l = 0; l < LOAD_REPEAT_NUM; l++)
	{
		int blkdItNone = 0, blkdItConv = 0, blkdItProp = 0, blkdItConvIlp = 0, blkdItPropIlp = 0;
		int togOpItConv = 0, togOpItProp = 0;
		int realOpItConv = 0, realOpItProp = 0;
		int rerouteOpItProp = 0;
		double stdItNoneArr[ITERATION] = {}, stdItConvArr[ITERATION] = {}, stdItPropArr[ITERATION] = {}, stdItConvIlpArr[ITERATION] = {}, stdItPropIlpArr[ITERATION] = {};
		if(l != 0){
			load = load + LOAD_GAIN;
		}
		ofs_result_txt 	<< endl << "load= "<< load << endl;
		cout 			<< endl << "load= "<< load << endl;
		initialize();
		reInitialize();
		if (readInput(argc, &argv[0], load) == 1){
		    cout << "[error] cannot read input" << endl;
		    return 1;
		}
		
		for (int it = 0; it < ITERATION; it++)
		{
			ofs_result_txt 	<< endl << "Iteration= "<< it << endl;
			cout 			<< endl << "Iteration= "<< it << endl;
			seed1 = (it+1) * seed1;
			seed2 = (it+2) * seed1;

			genDemands(load);

			reInitialize();
			initializeEvent();

			for(algoCall = 0; algoCall <= 2; algoCall++)
			{
				for (int i = 0; i < REQ_NUM; i++)
				{
					eventQueue.push(endEvent[i]);
					eventQueue.push(startEvent[i]);
				}
				for (int c = 0; c < defragCount; c++){
					eventQueue.push(defragEvent[c]);
				}

				time_slot_now = 0; 
				int lp = 0;
				clock_t start, end;
				start = clock();
				
				while(!eventQueue.empty())
				{
					int b = 0;
					while(!deleteQueue.empty())
					{
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
						
					nowEvent = eventQueue.top();
	 				eventQueue.pop();
	 				time_slot_now = nowEvent.time_slot_now;
	 				if(nowEvent.type == 1){
	 					time_slot_now = nowEvent.time_slot_now;
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
	 				if(nowEvent.type == 0)
	 				{
	 					if(lp_size[nowEvent.lpNum])
	 					{
	 						last_lp = nowEvent.lpNum;
	 						try {
								b = firstFit1_1(nowEvent.lpNum, algoCall);
							}
							catch(const char* err){
								cout << "ERR:パス割り当て中における" << err << endl;
								return 1;
							}
							if(!b)
							{
								try {
									retuneBp(load);
								}
								catch(const char* err){
									cout << "ERR:ブロッキング中における" << err << endl;
									return 1;
								}
							}
							int s= source[nowEvent.lpNum], d= dest[nowEvent.lpNum], sort_val1 = 0, sort_val2 = 0;
							if(b) addToList(nowEvent.lpNum, sort_val1, sort_val2);
							if (nowEvent.lpNum == REQ_NUM-1)
							{
								while(!eventQueue.empty())
								{
									eventQueue.pop();
								}
							}
	 					}
	 				}
					
					if(nowEvent.type == 2){
						try{
							retuneBp(load);
						}
						catch(const char* err)
						{
							cout << "ERR:デフラグ中における" << err << endl;
							return 1;
						}						
					}
				}
				
				end = clock();
   				cout << "Simulation time_slot_now = " << (double)(end - start) / CLOCKS_PER_SEC << " sec" <<endl;
   				ofs_result_txt << "Simulation time_slot_now = " << (double)(end - start) / CLOCKS_PER_SEC << " sec" <<endl;

				if(algoCall==0){
					cout 			<< "Blocked request: First-fit w/o defragment                " << blocked << endl;
					ofs_result_txt 	<< "Blocked request: First-fit w/o defragment                " << blocked << endl << endl;
					blkdItNone 			+= blocked;
					stdItNoneArr[it] 	= blocked;
					reInitialize();
				}
				if(algoCall==1){
					cout 			<< "Blocked request: First-fit w/ stat defragment            " << blocked << endl;
					ofs_result_txt 	<< "Blocked request: First-fit w/ stat defragment            " << blocked << endl << endl;
					cout 			<< "Number of switching operations:                          " << togOp << endl;
					ofs_result_txt	<< "Number of switching operations:                          " << togOp << endl;
					cout 			<< "Number of reallocating operations:                       " << realOp << endl;
					ofs_result_txt 	<< "Number of reallocating operations:                       " << realOp << endl << endl;
					blkdItConv 			+= blocked;
					stdItConvArr[it] 	= blocked;
					togOpItConv  		+= togOp;
					realOpItConv 		+= realOp;
					reInitialize();
				}
				if(algoCall==2){
					cout 			<< "Blocked request: First-fit w/ stat rerouting defragment  " << blocked << endl;
					ofs_result_txt 	<< "Blocked request: First-fit w/ stat rerouting defragment  " << blocked << endl << endl;
					cout 			<< "Number of switching operations:                          " << togOp << endl;
					ofs_result_txt 	<< "Number of switching operations:                          " << togOp << endl;
					cout 			<< "Number of reallocating operations:                       " << realOp << endl;
					ofs_result_txt 	<< "Number of reallocating operations:                       " << realOp << endl << endl;
					cout 			<< "Number of rerouting operations:                          " << rerouteOp << endl;
					ofs_result_txt 	<< "Number of rerouting operations:                          " << rerouteOp << endl << endl;
					blkdItProp 			+= blocked;
					stdItPropArr[it] 	= blocked;
					togOpItProp  		+= togOp;
					realOpItProp 		+= realOp;
					rerouteOpItProp 	+= rerouteOp;
					reInitialize();
				}
				if(algoCall==3){
					cout 			<< "Blocked request: First-fit w/ stat ilp defragment        " << blocked << endl;
					ofs_result_txt 	<< "Blocked request: First-fit w/ stat ilp defragment        " << blocked << endl << endl;
					cout 			<< "Number of switching operations:                          " << togOp << endl;
					ofs_result_txt 	<< "Number of switching operations:                          " << togOp << endl;
					cout 			<< "Number of reallocating operations:                       " << realOp << endl;
					ofs_result_txt 	<< "Number of reallocating operations:                       " << realOp << endl << endl;
					blkdItConvIlp 		+= blocked;
					stdItConvIlpArr[it] = blocked;
					reInitialize();
				}
				if(algoCall==4){
					cout 			<< "Blocked request: First-fit w/ stat ilp reroute defragment" << blocked << endl;
					ofs_result_txt 	<< "Blocked request: First-fit w/ stat ilp reroute defragment" << blocked << endl << endl;
					cout 			<< "Number of switching operations:                          " << togOp << endl;
					ofs_result_txt 	<< "Number of switching operations:                          " << togOp << endl;
					cout 			<< "Number of reallocating operations:                       " << realOp << endl;
					ofs_result_txt 	<< "Number of reallocating operations:                       " << realOp << endl << endl;
					blkdItPropIlp += blocked;
					stdItPropIlpArr[it] = blocked;
					reInitialize();
				}
			}
		}
		ofs_result_csv << load << ",";
		operation_num_csv << load << ",";


		double stdItNone, stdItConv, stdItProp, stdItConvIlp, stdItPropIlp;
		double stdItNoneDiff, stdItConvDiff, stdItPropDiff, stdItConvIlpDiff, stdItPropIlpDiff;
		stdItNone = standard(stdItNoneArr, ITERATION);
		stdItConv = standard(stdItConvArr, ITERATION);
		stdItProp = standard(stdItPropArr, ITERATION);
		stdItConvIlp = standard(stdItConvIlpArr, ITERATION);
		stdItPropIlp = standard(stdItPropIlpArr, ITERATION);
		stdItNoneDiff = blkdItNone*0.05 - stdItNone*1.96;
		stdItConvDiff = blkdItConv*0.05 - stdItConv*1.96;
		stdItPropDiff = blkdItProp*0.05 - stdItProp*1.96;
		stdItConvIlpDiff = blkdItConvIlp*0.05 - stdItConvIlp*1.96;
		stdItPropIlpDiff = blkdItPropIlp*0.05 - stdItPropIlp*1.96;

		cout << "Av blocked request blkdItNone:         " << blkdItNone/ITERATION << endl;
		ofs_result_txt << "Av blocked request blkdItNone:         " << blkdItNone/ITERATION << endl;
		ofs_result_csv << blkdItNone/ITERATION << ",";
		cout << "Confidence interval blkdItNone:        " << stdItNone << endl;
		ofs_result_txt << "Confidence interval blkdItNone:        " << stdItNone << endl;
		cout << "Av blocked request blkdItConv:         " << blkdItConv/ITERATION << endl;
		ofs_result_txt << "Av blocked request blkdItConv:         " << blkdItConv/ITERATION << endl;
		ofs_result_csv << blkdItConv/ITERATION << ",";
		cout << "Confidence interval blkdItConv:        " << stdItConv << endl;
		ofs_result_txt << "Confidence interval blkdItConv:        " << stdItConv << endl;
		cout << "Av blocked request blkdItProp:         " << blkdItProp/ITERATION << endl;
		ofs_result_txt << "Av blocked request blkdItProp:         " << blkdItProp/ITERATION << endl;
		ofs_result_csv << blkdItProp/ITERATION << ",";
		cout << "Confidence interval blkdItProp:        " << stdItProp << endl;
		ofs_result_txt << "Confidence interval blkdItProp:        " << stdItProp << endl;
		cout << "Av blocked request blkdItConvIlp:      " << blkdItConvIlp/ITERATION << endl;
		ofs_result_txt << "Av blocked request blkdItConvIlp:      " << blkdItConvIlp/ITERATION << endl;
		ofs_result_csv << blkdItConvIlp/ITERATION << ",";
		cout << "Confidence interval blkdItConvIlp:     " << stdItConvIlp << endl;
		ofs_result_txt << "Confidence interval blkdItConvIlp:     " << stdItConvIlp << endl;
		cout << "Av blocked request blkdItPropIlp:      " << blkdItPropIlp/ITERATION << endl;
		ofs_result_txt << "Av blocked request blkdItPropIlp:      " << blkdItPropIlp/ITERATION << endl;
		ofs_result_csv << blkdItPropIlp/ITERATION << endl;
		cout << "Confidence interval blkdItPropIlp:     " << stdItPropIlp << endl;
		ofs_result_txt << "Confidence interval blkdItPropIlp:     " << stdItPropIlp << endl;
		cout << "Av toggle operations togOpItConv:      " << togOpItConv/ITERATION << endl;
		ofs_result_txt << "Av toggle operations togOpItConv:      " << togOpItConv/ITERATION << endl;
		operation_num_csv << togOpItConv/ITERATION << ",";
		cout << "Av move operations realOpItConv:       " << realOpItConv/ITERATION << endl;
		ofs_result_txt << "Av move operations realOpItConv:       " << realOpItConv/ITERATION << endl;
		operation_num_csv << realOpItConv/ITERATION << ",";
		cout << "Av toggle operations togOpItProp:      " << togOpItProp/ITERATION << endl;
		ofs_result_txt << "Av toggle operations togOpItProp:      " << togOpItProp/ITERATION << endl;
		operation_num_csv << togOpItProp/ITERATION << ",";
		cout << "Av move operations realOpItProp:       " << realOpItProp/ITERATION << endl;
		ofs_result_txt << "Av move operations realOpItProp:       " << realOpItProp/ITERATION << endl;
		operation_num_csv << realOpItProp/ITERATION << ",";
		cout << "Av reroute operations rerouteOpItProp: " << rerouteOpItProp/ITERATION << endl << endl;
		ofs_result_txt << "Av reroute operations rerouteOpItProp: " << rerouteOpItProp/ITERATION << endl << endl;
		operation_num_csv << rerouteOpItProp/ITERATION << endl;
	}
	ofs_result_txt.close();

	lpNode *cur;
	cur = realList;
	while (cur != NULL) 
	{
		deleteLP(cur->x, 0);
		cur = cur->next;
	}

	cout << endl;
	cout << "time_slot_now: " << time_slot_now << " Seed1: " << seed1 << ", Seed2: " << seed2 << endl<< endl;

	return 0;
}