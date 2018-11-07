#define INF 999999999
#define ITERATION 3
#define REQ_NUM 100
#define NODE_NUM 11

// directional links (NODE_NUM,LINK_NUM) = (11, 28), (5, 12), (14,44), (11, 52), (14, 46), (25, 84)
#define LINK_NUM 52

#define CAPASITY 400
#define START_LOAD 350

// average holding time
#define HOLDING_TIME 10

// max demand size
#define req_Max 16

// fixme:smallness but it is not used
#define K 1000

#define DEFRAG_TIME 0.1

// maximum total defrag time
#define  T_temp 100

// Retuning period
#define R_int 0.2

// sort type
// 0 if not used, 1 by length, 2 by size, 3 by blocks
// #define Sorting 1

// maxSteo for ILP
#define maxStep 3

// defrag speed but it is not used
#define spfact 1

// max hop num
#define LIMIT_HOP_NUM 99999