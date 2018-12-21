# defragmentation simulator

## Requirements

- gcc 4.2.1
- IBM(R) ILOG(R) CPLEX(R) Interactive Optimizer 12.7.1.0
- Python 3.5.2

## Usage

```
# download simulator
git clone https://github.com/adshidtadka/defragmentation

# create result folder
cd defragmentation
mkdir result

# move to apprecation folder
cd application

# put parameter
vim params.h
## choose topology 
## from (NODE_NUM,LINK_NUM) = (11, 28), (5, 12), (14,44), (11, 52), (14, 46), (25, 84)
## and put them respectively

# compile
make

# run simulator
./main.out [5node]
## choose topology
## from (5node, abi, euro, nsf, jap)
```


## Notes

|  topology  |  5node  |  abi  |  nsf  |  euro  |  jap  |
| ---- | ---- | ---- | ---- | ---- | ---- |
|  NODE_NUM  |  5   |  11  |  14  |  11  |  25  |
|  LINK_NUM  |  12  |  28  |  44  |  52  |  84  |




```
