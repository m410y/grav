#ifndef PARTICLE_H
#define PARTICLE_H

#include <metric_field.h>

typedef struct particle Particle;
struct particle {
    Rank1 pos;
    Rank1 vel;
    double time;
};

void move_particle(Particle *particle, Metric_field *field, double ds);

#endif
