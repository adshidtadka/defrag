#define INF 999999999
#define T_end 1245201000

// iteration
#define IT 3

// demand
#define M 10000

// node num
#define N 11

// link num 
// directional links (N,L) = (11, 28), (5, 12), (14,44), (11, 52), (14, 46), (25, 84)
#define L 52

// capacity per link
#define S 400

// start traffic load
#define A1 350

// average holding time
#define H  10

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