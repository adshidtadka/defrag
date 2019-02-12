struct Constant {
    static const int INF                            = 999999999;
    static const int ITERATION                      = 3;
    static const int REQ_NUM                        = 100;
    static const int CAPACITY                       = 400;
    static const int LOAD_START                     = 50;
    static const int LOAD_GAIN                      = 10;
    static const int LOAD_REPEAT_NUM		        = 10;
    static const int HOLDING_TIME 			        = 10;
    static const int REQ_SIZE_MAX 			        = 16;
    constexpr static const double PROCESSING_TIME 	= 0.1;
    constexpr static const double DEFRAG_INTERVAL 	= 0.2;
    static const int DEFRAG_TOTAL_TIME_MAX          = 100;
    static const int MAX_STEP 				        = 3;
    static const int ADDITIONAL_HOP 		        = 0;
    static const int ADDITIONAL_HOP_CONGESTED       = 10;
    static const int RATE_CONGESTED                 = 20;
    static const int SEED_1				            = 1448601515;
    static const int SEED_2				            = 125;
    static const int ALGO_START			            = 0;
    static const int ALGO_FINISH			        = 2;
    static const int NODE_NUM 				        = 25;
    static const int LINK_NUM 				        = 84;
    // (NODE_NUM,LINK_NUM) = (11, 28), (5, 12), (14,44), (11, 52), (14, 46), (25, 84)
};
