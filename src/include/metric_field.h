#ifndef METRIC_FIELD_H
#define METRIC_FIELD_H

#include <tensor_deifnes.h>

#include <stddef.h>

typedef struct face Face;
struct face {
    Face *near[3];
    Rank1 *vert[3];
    Rank2 metric;
    Rank3 d_metric;
};

typedef struct metric_field Metric_field;
struct metric_field {
    Face *f;
    Rank1 *vertices;
};

typedef Rank2 (*metric_func)(Rank1);

int field_alloc(Metric_field *field, size_t Nx, size_t Ny);
void field_delete(Metric_field *field);
void field_func_init(Metric_field *field, metric_func func);
Rank2 metric_at_point(Metric_field *field, Rank1 P);
Rank3 christoffel_at_point(Metric_field *field, Rank1 P);

#endif
