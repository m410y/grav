#ifndef UNIT_TESTS

#include <particle.h>

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

Rank2 metric(Rank1 point)
{
    const double x = point.x - 0.5;
    const double y = point.y - 0.5;
    const double z = 10 * exp(-(x * x + y * y) / 2);
    const Rank2 metric = {
        .x.x = 1 + x * x * z,
        .x.y = x * y * z,
        .y.x = y * x * z,
        .y.y = 1 + y * y * z,
    };
    return metric;
}

void particle_init(Particle *particle, FILE *const file)
{
    fscanf(file, "%lf %lf", &particle->pos.x, &particle->pos.y);
    fscanf(file, "%lf %lf", &particle->vel.x, &particle->vel.y);

    particle->time = 0;
}

void particle_log(Particle *particle, FILE *const file)
{
    fprintf(file, "%lf\t%lf\t%lf\t%lf\t%lf\n", particle->time,
            particle->pos.x, particle->pos.y,
            particle->vel.x, particle->vel.y);
}

int main()
{
    int i, max_i;
    size_t Nx, Ny;
    double step;
    FILE *in, *out;
    Metric_field field;
    Particle tester;

    in = fopen("resources/input.txt", "r");
    fscanf(in, "%d", &max_i);
    fscanf(in, "%lf", &step);
    fscanf(in, "%lu%lu", &Nx, &Ny);

    if (field_alloc(&field, Nx, Ny) < 0)
    {
        fclose(in);
        return -1;
    }
    field_func_init(&field, metric);

    particle_init(&tester, in);
    fclose(in);

    out = fopen("resources/output.txt", "w");
    particle_log(&tester, out);
    for (i = 0; i < max_i; i++)
    {
        move_particle(&tester, &field, step);
        particle_log(&tester, out);
    }

    fclose(out);
    field_delete(&field);
    return 0;
}

#endif
