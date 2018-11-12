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
		ofs_result_txt 	<< endl << "load = "<< load << endl;
		cout 			<< endl << "load = "<< load << endl;
		initialize();
		if (readInput(argc, &argv[0], load) == 1){
		    cout << "[error] cannot read input" << endl;
		    return 1;
		}
		
		for (int it = 0; it < ITERATION; it++)
		{
			ofs_result_txt 	<< endl << "Iteration = "<< it << endl;
			cout 			<< endl << "Iteration = "<< it << endl;
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
					    	removeLP1_1(nowEvent.lpNum);
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
					    	removeLP1_1(nowEvent.lpNum);
					    }
					    catch(const char* err) {
					    	cout << "ERR:パス切断中における" << err << endl;
					    	return 1;
					    }
						delFromList(1, nowEvent.lpNum);
						delFromList(2, nowEvent.lpNum);
	 				}
	 				if(nowEvent.type == 0)
	 				{
	 					if(lp_size[nowEvent.lpNum])
	 					{
	 						last_lp = nowEvent.lpNum;
	 						try {
								b = firstFit1_1(nowEvent.lpNum);
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
							int sort_val1 = 0, sort_val2 = 0;
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
   				cout 				<< "[algo"<< algoCall << "] it takes " << (double)(end - start) / CLOCKS_PER_SEC << " sec" <<endl;
   				ofs_result_txt 		<< "[algo"<< algoCall << "] it takes " << (double)(end - start) / CLOCKS_PER_SEC << " sec" <<endl;

				if(algoCall==0){
					cout 			<< "[algo0] blocked:		" << blocked << endl;
					ofs_result_txt 	<< "[algo0] blocked:		" << blocked << endl << endl;
					blkdItNone 			+= blocked;
					stdItNoneArr[it] 	= blocked;
					reInitialize();
				}
				if(algoCall==1){
					cout 			<< "[algo1] blocked:		" << blocked << endl;
					ofs_result_txt 	<< "[algo1] blocked:		" << blocked << endl << endl;
					cout 			<< "[algo1] toggle:		" << togOp << endl;
					ofs_result_txt	<< "[algo1] toggle:		" << togOp << endl;
					cout 			<< "[algo1] reallocate:	" << realOp << endl;
					ofs_result_txt 	<< "[algo1] reallocate:	" << realOp << endl << endl;
					blkdItConv 			+= blocked;
					stdItConvArr[it] 	= blocked;
					togOpItConv  		+= togOp;
					realOpItConv 		+= realOp;
					reInitialize();
				}
				if(algoCall==2){
					cout 			<< "[algo2] blocked:		" << blocked << endl;
					ofs_result_txt 	<< "[algo2] blocked:		" << blocked << endl << endl;
					cout 			<< "[algo2] toggle:		" << togOp << endl;
					ofs_result_txt 	<< "[algo2] toggle:		" << togOp << endl;
					cout 			<< "[algo2] reallocate:	" << realOp << endl;
					ofs_result_txt 	<< "[algo2] reallocate:	" << realOp << endl << endl;
					cout 			<< "[algo2] reroute:	" << rerouteOp << endl;
					ofs_result_txt 	<< "[algo2] reroute:	" << rerouteOp << endl << endl;
					blkdItProp 			+= blocked;
					stdItPropArr[it] 	= blocked;
					togOpItProp  		+= togOp;
					realOpItProp 		+= realOp;
					rerouteOpItProp 	+= rerouteOp;
					reInitialize();
				}
				if(algoCall==3){
					cout 			<< "[algo3] blocked:		" << blocked << endl;
					ofs_result_txt 	<< "[algo3] blocked:		" << blocked << endl << endl;
					cout 			<< "[algo3] toggle:		" << togOp << endl;
					ofs_result_txt 	<< "[algo3] toggle:		" << togOp << endl;
					cout 			<< "[algo3] reallocate:	" << realOp << endl;
					ofs_result_txt 	<< "[algo3] reallocate:	" << realOp << endl << endl;
					blkdItConvIlp 		+= blocked;
					stdItConvIlpArr[it] = blocked;
					reInitialize();
				}
				if(algoCall==4){
					cout 			<< "[algo4] blocked:		" << blocked << endl;
					ofs_result_txt 	<< "[algo4] blocked:		" << blocked << endl << endl;
					cout 			<< "[algo4] toggle:		" << togOp << endl;
					ofs_result_txt 	<< "[algo4] toggle:		" << togOp << endl;
					cout 			<< "[algo4] reallocate:	" << realOp << endl;
					ofs_result_txt 	<< "[algo4] reallocate:	" << realOp << endl << endl;
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

		cout				<< endl << "Average results" << endl;
		ofs_result_txt		<< endl << "Average results" << endl;
		cout 				<< "[algo0] blocked:		" << blkdItNone/ITERATION << endl;
		ofs_result_txt 		<< "[algo0] blocked:		" << blkdItNone/ITERATION << endl;
		ofs_result_csv 		<< blkdItNone/ITERATION << ",";
		cout 				<< "[algo0] confidence:	" << stdItNone << endl << endl;
		ofs_result_txt 		<< "[algo0] confidence:	" << stdItNone << endl << endl;
		cout 				<< "[algo1] blocked:		" << blkdItConv/ITERATION << endl;
		ofs_result_txt 		<< "[algo1] blocked:		" << blkdItConv/ITERATION << endl;
		ofs_result_csv 		<< blkdItConv/ITERATION << ",";
		cout 				<< "[algo1] confidence:	" << stdItConv << endl << endl;
		ofs_result_txt 		<< "[algo1] confidence:	" << stdItConv << endl << endl;
		cout 				<< "[algo2] blocked:		" << blkdItProp/ITERATION << endl;
		ofs_result_txt 		<< "[algo2] blocked:		" << blkdItProp/ITERATION << endl;
		ofs_result_csv 		<< blkdItProp/ITERATION << ",";
		cout 				<< "[algo2] confidence:	" << stdItProp << endl << endl;
		ofs_result_txt 		<< "[algo2] confidence:	" << stdItProp << endl << endl;
		cout 				<< "[algo3] blocked:		" << blkdItConvIlp/ITERATION << endl;
		ofs_result_txt 		<< "[algo3] blocked:		" << blkdItConvIlp/ITERATION << endl;
		ofs_result_csv 		<< blkdItConvIlp/ITERATION << ",";
		cout 				<< "[algo3] confidence:	" << stdItConvIlp << endl << endl;
		ofs_result_txt 		<< "[algo3] confidence:	" << stdItConvIlp << endl << endl;
		cout 				<< "[algo4] blocked:		" << blkdItPropIlp/ITERATION << endl;
		ofs_result_txt 		<< "[algo4] blocked:		" << blkdItPropIlp/ITERATION << endl;
		ofs_result_csv 		<< blkdItPropIlp/ITERATION << endl;
		cout 				<< "[algo4] confidence:	" << stdItPropIlp << endl << endl;
		ofs_result_txt 		<< "[algo4] confidence:	" << stdItPropIlp << endl << endl;
		cout 				<< "[algo1] toggle:		" << togOpItConv/ITERATION << endl;
		ofs_result_txt 		<< "[algo1] toggle:		" << togOpItConv/ITERATION << endl;
		operation_num_csv 	<< togOpItConv/ITERATION << ",";
		cout 				<< "[algo1] reallocate:	" << realOpItConv/ITERATION << endl << endl;
		ofs_result_txt 		<< "[algo1] reallocate:	" << realOpItConv/ITERATION << endl << endl;
		operation_num_csv 	<< realOpItConv/ITERATION << ",";
		cout 				<< "[algo2] toggle:		" << togOpItProp/ITERATION << endl;
		ofs_result_txt 		<< "[algo2] toggle:		" << togOpItProp/ITERATION << endl;
		operation_num_csv 	<< togOpItProp/ITERATION << ",";
		cout 				<< "[algo2] reallocate:	" << realOpItProp/ITERATION << endl;
		ofs_result_txt 		<< "[algo2] reallocate:	" << realOpItProp/ITERATION << endl;
		operation_num_csv 	<< realOpItProp/ITERATION << ",";
		cout 				<< "[algo2] reroute:	" << rerouteOpItProp/ITERATION << endl << endl;
		ofs_result_txt 		<< "[algo2] reroute:	" << rerouteOpItProp/ITERATION << endl << endl;
		operation_num_csv 	<< rerouteOpItProp/ITERATION << endl;
	}
	ofs_result_txt.close();

	cout << endl;
	cout << "time_slot_now: " << time_slot_now << " Seed1: " << seed1 << ", Seed2: " << seed2 << endl<< endl;

	return 0;
}