#ifndef TENSOR_DEFINES_H
#define TENSOR_DEFINES_H

typedef double Rank0;

typedef struct rank1 {
    Rank0 x, y;
} Rank1;

typedef struct rank2 {
    Rank1 x, y;
} Rank2;

typedef struct rank3 {
    Rank2 x, y;
} Rank3;

typedef struct rank4 {
    Rank3 x, y;
} Rank4;

#define PTR_CAST(var) ((double *)&(var))

#endif
