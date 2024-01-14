#ifndef METRIC_FIELD_H
#define METRIC_FIELD_H

#include <tensor_deifnes.h>

typedef struct metric_field Metric_field;
struct metric_field {
    void *f;
    Rank1 *vertices;
};

typedef Rank2 (*metric_func)(Rank1);

int field_alloc(Metric_field *field, int Nx, int Ny);
void field_delete(Metric_field *field);
void field_func_init(Metric_field *field, metric_func func);
Rank2 field_metric_at_point(Metric_field *field, Rank1 P);
Rank3 field_christoffel_at_point(Metric_field *field, Rank1 P);
Rank4 field_riemann_at_point(Metric_field *field, Rank1 P);
Rank2 field_ricci_at_point(Metric_field *field, Rank1 P);
Rank0 field_curavture_at_point(Metric_field *field, Rank1 P);

#endif
