#include "params.h"
#include "func.h"

int main(int argc, char *argv[]) {
    ofstream ofs_result_txt;
    ofs_result_txt.open("./../result/result.txt", ios::out);
    ofstream ofs_result_csv;
    ofs_result_csv.open("./../result/result.csv", ios::out);
    ofstream ofs_result_operation_csv;
    ofs_result_operation_csv.open("./../result/result_operation.csv", ios::out);
    ofstream ofs_slot_csv;
    ofs_slot_csv.open("./../result/slot.csv", ios::out);

    ofs_result_txt << "-------- Simulation start --------" << endl;
    ofs_result_txt << "Constant::NODE_NUM 		= " << Constant::NODE_NUM << endl;
    ofs_result_txt << "Constant::LINK_NUM 		= " << Constant::LINK_NUM << endl;
    ofs_result_txt << "Constant::CAPACITY 		= " << Constant::CAPACITY << endl;
    ofs_result_txt << "Constant::REQ_NUM 		= " << Constant::REQ_NUM << endl;
    ofs_result_txt << "Constant::REQ_SIZE_MAX 		= " << Constant::REQ_SIZE_MAX << endl;
    ofs_result_txt << "Constant::DEFRAG_INTERVAL 	= " << Constant::DEFRAG_INTERVAL << endl;
    ofs_result_txt << "Constant::DEFRAG_TOTAL_TIME_MAX 	= " << Constant::DEFRAG_TOTAL_TIME_MAX << endl;
    ofs_result_txt << "Constant::PROCESSING_TIME 	= " << Constant::PROCESSING_TIME << endl;
    ofs_result_txt << "Constant::ADDITIONAL_HOP 		= " << Constant::ADDITIONAL_HOP << endl;
    cout << "-------- Simulation start --------" << endl;
    cout << "Constant::NODE_NUM 		= " << Constant::NODE_NUM << endl;
    cout << "Constant::LINK_NUM 		= " << Constant::LINK_NUM << endl;
    cout << "Constant::CAPACITY 		= " << Constant::CAPACITY << endl;
    cout << "Constant::REQ_NUM 		= " << Constant::REQ_NUM << endl;
    cout << "Constant::REQ_SIZE_MAX 		= " << Constant::REQ_SIZE_MAX << endl;
    cout << "Constant::DEFRAG_INTERVAL 	= " << Constant::DEFRAG_INTERVAL << endl;
    cout << "Constant::DEFRAG_TOTAL_TIME_MAX	= " << Constant::DEFRAG_TOTAL_TIME_MAX << endl;
    cout << "Constant::PROCESSING_TIME 	= " << Constant::PROCESSING_TIME << endl;
    cout << "Constant::ADDITIONAL_HOP 	= " << Constant::ADDITIONAL_HOP << endl;

    int load = Constant::LOAD_START;
    for (int l = 0; l < Constant::LOAD_REPEAT_NUM; l++) {
        int blkdItNone = 0, blkdItConv = 0, blkdItProp = 0, blkdItConvIlp = 0, blkdItPropIlp = 0;
        int togOpItConv = 0, togOpItProp = 0;
        int realOpItConv = 0, realOpItProp = 0;
        int rerouteOpItProp = 0;
        double stdItNoneArr[Constant::ITERATION] = {}, stdItConvArr[Constant::ITERATION] = {}, stdItPropArr[Constant::ITERATION] = {}, stdItConvIlpArr[Constant::ITERATION] = {}, stdItPropIlpArr[Constant::ITERATION] = {};
        int slotNone[Constant::LINK_NUM][Constant::ITERATION], slotConv[Constant::LINK_NUM][Constant::ITERATION], slotProp[Constant::LINK_NUM][Constant::ITERATION], slotConvIlp[Constant::LINK_NUM][Constant::ITERATION], slotPropIlp[Constant::LINK_NUM][Constant::ITERATION];
        if (l != 0) {
            load = load + Constant::LOAD_GAIN;
        }
        ofs_result_txt << endl << "---- load = " << load << " ----" << endl;
        cout << endl << "---- load = " << load << " ----" << endl;
        initialize();
        if (readInput(argc, &argv[0], load) == 1) {
            cout << "[error] cannot read input" << endl;
            return 1;
        }

        congestedLink();

        for (int it = 0; it < Constant::ITERATION; it++) {
            ofs_result_txt << endl << "Iteration = " << it << endl;
            cout << endl << "Iteration = " << it << endl;
            seed1 = (it + 1) * seed1;
            seed2 = (it + 2) * seed1;

            genDemands(load);
            reInitialize();
            initializeEvent();

            for (algoCall = Constant::ALGO_START; algoCall <= Constant::ALGO_FINISH; algoCall++) {
                for (int i = 0; i < Constant::REQ_NUM; i++) {
                    eventQueue.push(endEvent[i]);
                    eventQueue.push(startEvent[i]);
                }
                for (int c = 0; c < defragCount; c++) {
                    eventQueue.push(defragEvent[c]);
                }

                vector <vector<int>> allocatedSlotNum(Constant::LINK_NUM);

                time_slot_now = 0;
                clock_t start, end;
                start = clock();

                while (!eventQueue.empty()) {
                    int b = 0;
                    while (!deleteQueue.empty()) {
                        nowEvent = deleteQueue.top();
                        deleteQueue.pop();
                        try {
                            removeLP1_1(nowEvent.lpNum);
                        }
                        catch (const char *err) {
                            cout << "ERR:パス切断中における" << err << endl;
                            return 1;
                        }
                        delFromList(1, nowEvent.lpNum);
                        delFromList(2, nowEvent.lpNum);
                    }

                    nowEvent = eventQueue.top();
                    eventQueue.pop();
                    time_slot_now = nowEvent.time;

                    if (nowEvent.type == 1) {
                        time_slot_now = nowEvent.time;
                        try {
                            removeLP1_1(nowEvent.lpNum);
                        }
                        catch (const char *err) {
                            cout << "ERR:パス切断中における" << err << endl;
                            return 1;
                        }
                        delFromList(1, nowEvent.lpNum);
                        delFromList(2, nowEvent.lpNum);
                    }

                    if (nowEvent.type == 0) {
                        if (lp_size[nowEvent.lpNum]) {
                            try {
                                b = setupPath(nowEvent.lpNum);
                            }
                            catch (const char *err) {
                                cout << "ERR:パス割り当て中における" << err << endl;
                                return 1;
                            }
                            if (!b) {
                                try {
                                    startDefrag(load);
                                }
                                catch (const char *err) {
                                    cout << "ERR:ブロッキング中における" << err << endl;
                                    return 1;
                                }
                            }
                            int sort_val1 = 0, sort_val2 = 0;
                            if (b) addToList(nowEvent.lpNum, sort_val1, sort_val2);
                            if (nowEvent.lpNum == Constant::REQ_NUM - 1) {
                                while (!eventQueue.empty()) {
                                    eventQueue.pop();
                                }
                            }
                        }
                    }

                    if (nowEvent.type == 2) {
                        try {
                            startDefrag(load);
                        }
                        catch (const char *err) {
                            cout << "ERR:デフラグ中における" << err << endl;
                            return 1;
                        }
                    }

                    // count allocated slot number
                    for (int i = 0; i < Constant::LINK_NUM; ++i) {
                        int counter = 0;
                        for (int j = 0; j < Constant::CAPACITY; ++j) {
                            if (spec[j][i]) {
                                counter++;
                            }
                        }
                        allocatedSlotNum[i].push_back(counter);
                    }
                }

                end = clock();
                cout << "algoCall = " << algoCall << ", it takes " << (double) (end - start) / CLOCKS_PER_SEC << " sec"
                     << endl;
                ofs_result_txt << "algoCall = " << algoCall << ", it takes " << (double) (end - start) / CLOCKS_PER_SEC
                               << " sec" << endl;

                if (algoCall == 0) {
                    cout << "[none] blocked:		" << blocked << endl;
                    ofs_result_txt << "[none] blocked:		" << blocked << endl << endl;
                    blkdItNone += blocked;
                    stdItNoneArr[it] = blocked;
                    for (int i = 0; i < Constant::LINK_NUM; ++i) {
                        slotNone[i][it] = accumulate(allocatedSlotNum[i].begin(), allocatedSlotNum[i].end(), 0) / allocatedSlotNum[i].size();
                    }
                    reInitialize();
                }
                if (algoCall == 1) {
                    cout << "[algo_conv] blocked:		" << blocked << endl;
                    ofs_result_txt << "[algo_conv] blocked:		" << blocked << endl << endl;
                    cout << "[algo_conv] toggle:		" << togOp << endl;
                    ofs_result_txt << "[algo_conv] toggle:		" << togOp << endl;
                    cout << "[algo_conv] reallocate:	" << realOp << endl;
                    ofs_result_txt << "[algo_conv] reallocate:	" << realOp << endl << endl;
                    blkdItConv += blocked;
                    stdItConvArr[it] = blocked;
                    togOpItConv += togOp;
                    realOpItConv += realOp;
                    for (int i = 0; i < Constant::LINK_NUM; ++i) {
                        slotConv[i][it] = accumulate(allocatedSlotNum[i].begin(), allocatedSlotNum[i].end(), 0) / allocatedSlotNum[i].size();
                    }
                    reInitialize();
                }
                if (algoCall == 2) {
                    cout << "[algo_prop] blocked:		" << blocked << endl;
                    ofs_result_txt << "[algo_prop] blocked:		" << blocked << endl << endl;
                    cout << "[algo_prop] toggle:		" << togOp << endl;
                    ofs_result_txt << "[algo_prop] toggle:		" << togOp << endl;
                    cout << "[algo_prop] reallocate:	" << realOp << endl;
                    ofs_result_txt << "[algo_prop] reallocate:	" << realOp << endl << endl;
                    cout << "[algo_prop] reroute:	" << rerouteOp << endl;
                    ofs_result_txt << "[algo_prop] reroute:	" << rerouteOp << endl << endl;
                    blkdItProp += blocked;
                    stdItPropArr[it] = blocked;
                    togOpItProp += togOp;
                    realOpItProp += realOp;
                    rerouteOpItProp += rerouteOp;
                    for (int i = 0; i < Constant::LINK_NUM; ++i) {
                        slotProp[i][it] = accumulate(allocatedSlotNum[i].begin(), allocatedSlotNum[i].end(), 0) / allocatedSlotNum[i].size();
                    }
                    reInitialize();
                }
                if (algoCall == 3) {
                    cout << "[ilp_conv] blocked:		" << blocked << endl;
                    ofs_result_txt << "[ilp_conv] blocked:		" << blocked << endl << endl;
                    cout << "[ilp_conv] toggle:		" << togOp << endl;
                    ofs_result_txt << "[ilp_conv] toggle:		" << togOp << endl;
                    cout << "[ilp_conv] reallocate:	" << realOp << endl;
                    ofs_result_txt << "[ilp_conv] reallocate:	" << realOp << endl << endl;
                    blkdItConvIlp += blocked;
                    stdItConvIlpArr[it] = blocked;
                    for (int i = 0; i < Constant::LINK_NUM; ++i) {
                        slotConvIlp[i][it] = accumulate(allocatedSlotNum[i].begin(), allocatedSlotNum[i].end(), 0) / allocatedSlotNum[i].size();
                    }
                    reInitialize();
                }
                if (algoCall == 4) {
                    cout << "[ilp_prop] blocked:		" << blocked << endl;
                    ofs_result_txt << "[ilp_prop] blocked:		" << blocked << endl << endl;
                    cout << "[ilp_prop] toggle:		" << togOp << endl;
                    ofs_result_txt << "[ilp_prop] toggle:		" << togOp << endl;
                    cout << "[ilp_prop] reallocate:	" << realOp << endl;
                    ofs_result_txt << "[ilp_prop] reallocate:	" << realOp << endl << endl;
                    blkdItPropIlp += blocked;
                    stdItPropIlpArr[it] = blocked;
                    for (int i = 0; i < Constant::LINK_NUM; ++i) {
                        slotPropIlp[i][it] = accumulate(allocatedSlotNum[i].begin(), allocatedSlotNum[i].end(), 0) / allocatedSlotNum[i].size();
                    }
                    reInitialize();
                }
            }
        }
        ofs_result_csv << load << ",";
        ofs_result_operation_csv << load << ",";

        double stdItNone, stdItConv, stdItProp, stdItConvIlp, stdItPropIlp;
        double stdItNoneDiff, stdItConvDiff, stdItPropDiff, stdItConvIlpDiff, stdItPropIlpDiff;
        stdItNone = standard(stdItNoneArr, Constant::ITERATION);
        stdItConv = standard(stdItConvArr, Constant::ITERATION);
        stdItProp = standard(stdItPropArr, Constant::ITERATION);
        stdItConvIlp = standard(stdItConvIlpArr, Constant::ITERATION);
        stdItPropIlp = standard(stdItPropIlpArr, Constant::ITERATION);
        stdItNoneDiff = blkdItNone * 0.05 - stdItNone * 1.96;
        stdItConvDiff = blkdItConv * 0.05 - stdItConv * 1.96;
        stdItPropDiff = blkdItProp * 0.05 - stdItProp * 1.96;
        stdItConvIlpDiff = blkdItConvIlp * 0.05 - stdItConvIlp * 1.96;
        stdItPropIlpDiff = blkdItPropIlp * 0.05 - stdItPropIlp * 1.96;

        cout << endl << "Average results" << endl;
        ofs_result_txt << endl << "Average results" << endl;
        cout << "[none] blocked:		" << blkdItNone / Constant::ITERATION << endl;
        ofs_result_txt << "[none] blocked:		" << blkdItNone / Constant::ITERATION << endl;
        ofs_result_csv << blkdItNone / Constant::ITERATION << ",";
        cout << "[none] confidence:	" << stdItNone << endl << endl;
        ofs_result_txt << "[none] confidence:	" << stdItNone << endl << endl;
        cout << "[algo_conv] blocked:		" << blkdItConv / Constant::ITERATION << endl;
        ofs_result_txt << "[algo_conv] blocked:		" << blkdItConv / Constant::ITERATION << endl;
        ofs_result_csv << blkdItConv / Constant::ITERATION << ",";
        cout << "[algo_conv] confidence:	" << stdItConv << endl << endl;
        ofs_result_txt << "[algo_conv] confidence:	" << stdItConv << endl << endl;
        cout << "[algo_prop] blocked:		" << blkdItProp / Constant::ITERATION << endl;
        ofs_result_txt << "[algo_prop] blocked:		" << blkdItProp / Constant::ITERATION << endl;
        ofs_result_csv << blkdItProp / Constant::ITERATION << ",";
        cout << "[algo_prop] confidence:	" << stdItProp << endl << endl;
        ofs_result_txt << "[algo_prop] confidence:	" << stdItProp << endl << endl;
        cout << "[ilp_conv] blocked:		" << blkdItConvIlp / Constant::ITERATION << endl;
        ofs_result_txt << "[ilp_conv] blocked:		" << blkdItConvIlp / Constant::ITERATION << endl;
        ofs_result_csv << blkdItConvIlp / Constant::ITERATION << ",";
        cout << "[ilp_conv] confidence:	" << stdItConvIlp << endl << endl;
        ofs_result_txt << "[ilp_conv] confidence:	" << stdItConvIlp << endl << endl;
        cout << "[ilp_prop] blocked:		" << blkdItPropIlp / Constant::ITERATION << endl;
        ofs_result_txt << "[ilp_prop] blocked:		" << blkdItPropIlp / Constant::ITERATION << endl;
        ofs_result_csv << blkdItPropIlp / Constant::ITERATION << endl;
        cout << "[ilp_prop] confidence:	" << stdItPropIlp << endl << endl;
        ofs_result_txt << "[ilp_prop] confidence:	" << stdItPropIlp << endl << endl;
        cout << "[algo_conv] toggle:		" << togOpItConv / Constant::ITERATION << endl;
        ofs_result_txt << "[algo_conv] toggle:		" << togOpItConv / Constant::ITERATION << endl;
        ofs_result_operation_csv << togOpItConv / Constant::ITERATION << ",";
        cout << "[algo_conv] reallocate:	" << realOpItConv / Constant::ITERATION << endl << endl;
        ofs_result_txt << "[algo_conv] reallocate:	" << realOpItConv / Constant::ITERATION << endl << endl;
        ofs_result_operation_csv << realOpItConv / Constant::ITERATION << ",";
        cout << "[algo_prop] toggle:		" << togOpItProp / Constant::ITERATION << endl;
        ofs_result_txt << "[algo_prop] toggle:		" << togOpItProp / Constant::ITERATION << endl;
        ofs_result_operation_csv << togOpItProp / Constant::ITERATION << ",";
        cout << "[algo_prop] reallocate:	" << realOpItProp / Constant::ITERATION << endl;
        ofs_result_txt << "[algo_prop] reallocate:	" << realOpItProp / Constant::ITERATION << endl;
        ofs_result_operation_csv << realOpItProp / Constant::ITERATION << ",";
        cout << "[algo_prop] reroute:	" << rerouteOpItProp / Constant::ITERATION << endl << endl;
        ofs_result_txt << "[algo_prop] reroute:	" << rerouteOpItProp / Constant::ITERATION << endl << endl;
        ofs_result_operation_csv << rerouteOpItProp / Constant::ITERATION << endl;
        for (int i = 0; i < Constant::LINK_NUM; ++i) {
            int slotNum = accumulate(slotNone[i], slotNone[i] + Constant::ITERATION, 0) / Constant::ITERATION;
            cout << "[none] allocated slot average number on link[" << i << "]: " << slotNum  << endl;
            ofs_slot_csv << slotNum << ",";
        }
        ofs_slot_csv << endl;
        for (int i = 0; i < Constant::LINK_NUM; ++i) {
            int slotNum = accumulate(slotConv[i], slotConv[i] + Constant::ITERATION, 0) / Constant::ITERATION;
            cout << "[conv] allocated slot average number on link[" << i << "]: " << slotNum  << endl;
            ofs_slot_csv << slotNum << ",";
        }
        ofs_slot_csv << endl;
        for (int i = 0; i < Constant::LINK_NUM; ++i) {
            int slotNum = accumulate(slotProp[i], slotProp[i] + Constant::ITERATION, 0) / Constant::ITERATION;
            cout << "[prop] allocated slot average number on link[" << i << "]: " << slotNum  << endl;
            ofs_slot_csv << slotNum << ",";
        }
        ofs_slot_csv << endl;
        for (int i = 0; i < Constant::LINK_NUM; ++i) {
            int slotNum = accumulate(slotConvIlp[i], slotConvIlp[i] + Constant::ITERATION, 0) / Constant::ITERATION;
            cout << "[conv_ilp] allocated slot average number on link[" << i << "]: " << slotNum  << endl;
            ofs_slot_csv << slotNum << ",";
        }
        ofs_slot_csv << endl;
        for (int i = 0; i < Constant::LINK_NUM; ++i) {
            int slotNum = accumulate(slotPropIlp[i], slotPropIlp[i] + Constant::ITERATION, 0) / Constant::ITERATION;
            cout << "[prop_ilp] allocated slot average number on link[" << i << "]: " << slotNum  << endl;
            ofs_slot_csv << slotNum << ",";
        }
        ofs_slot_csv << endl;
    }
    ofs_result_txt.close();

    cout << endl;
    cout << "time_slot_now: " << time_slot_now << " Seed1: " << seed1 << ", Seed2: " << seed2 << endl << endl;

    return 0;
}