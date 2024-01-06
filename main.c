#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct rank1 {
    double x, y;
} Rank1;

typedef struct rank2 {
    double xx, xy,
           yx, yy;
} Rank2;

typedef struct rank3 {
    double xxx, xxy, xyx, xyy,
           yxx, yxy, yyx, yyy;
} Rank3;

typedef struct rank4 {
    double xxxx, xxxy, xxyx, xxyy,
           xyxx, xyxy, xyyx, xyyy,

           yxxx, yxxy, yxyx, yxyy,
           yyxx, yyxy, yyyx, yyyy;
} Rank4;

typedef struct rank2_field {
    Rank2 *F;
    size_t Nx, Ny;
    double dx, dy;
} Rank2_field;

typedef struct particle_2d {
    Rank1 pos;
    Rank1 vel;
    double time;
} Particle_2d;

void rank2_field_alloc(Rank2_field *Field, size_t Nx, size_t Ny)
{
    if (Field->F)
        return;

    Field->F = malloc(Nx * Ny * sizeof(Rank2));
    if (!Field->F)
        return;

    Field->Nx = Nx;
    Field->Ny = Ny;
    Field->dx = 1.0 / Nx;
    Field->dy = 1.0 / Ny;
}

void rank2_field_init(Rank2_field *Field, Rank2 (*func)(Rank1))
{
    Rank1 point = {Field->dx / 2, Field->dy / 2};
    for (int i = 0; i < Field->Ny; i++)
    {
        for (int j = 0; j < Field->Nx; j++)
        {
            Field->F[i * Field->Nx + j] = func(point);
            point.x += Field->dx;
        }
        point.x = Field->dx / 2;
        point.y += Field->dy;
    }
}

void delete_field(Rank2_field *Field)
{
    free(Field->F);
}

Rank2 calc_metric(Rank1 point)
{
    Rank2 metric = {
        .xx = 1.0,
        .xy = 0.0,
        .yx = 0.0,
        .yy = sin(point.x * M_PI) * sin(point.x * M_PI),
    };
    return metric;
}

Rank2 rank2_inverse(Rank2 *M)
{
    double det = M->xx * M->yy - M->xy * M->yx;
    Rank2 res = {
        .xx = M->yy / det,
        .xy = -M->xy / det,
        .yx = -M->yx / det,
        .yy = M->xx / det,
    };

    return res;
}

Rank1 nearest_points(Rank2_field *Field, Rank1 *P, Rank2 *elems[4])
{
    int j = P->x / Field->dx - 0.5;
    int i = P->y / Field->dy - 0.5;
    Rank1 Pp = {
        .x = P->x * Field->Nx - j - 0.5,
        .y = P->y * Field->Ny - i - 0.5,
    };

    if (j < 0)
        j += Field->Nx;
    if (i < 0)
        i += Field->Ny;

    j %= Field->Nx;
    i %= Field->Ny;
    int jp = (j + 1) % Field->Nx;
    int ip = (i + 1) % Field->Ny;

    elems[0] = &Field->F[i * Field->Nx + j];
    elems[1] = &Field->F[i * Field->Nx + jp];
    elems[2] = &Field->F[ip * Field->Nx + j];
    elems[3] = &Field->F[ip * Field->Nx + jp];

    return Pp;
}

Rank2 value_at_point(Rank2_field *Field, Rank1 P)
{
    Rank2 *elems[4];
    P = nearest_points(Field, &P, elems);

    double coefs[4];
    coefs[3] = P.x * P.y;
    coefs[2] = P.y - coefs[3];
    coefs[0] = 1.0 - P.x - coefs[2];
    coefs[1] = P.x - coefs[3];

    Rank2 res;
    memset(&res, 0, sizeof(Rank2));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            ((double *)&res)[i] += coefs[j] *((double *)elems[j])[i];

    return res;
}

Rank3 diff_at_point(Rank2_field *Field, Rank1 P)
{
    Rank2 *elems[4];
    P = nearest_points(Field, &P, elems);

    double coefs_x[4];
    coefs_x[3] = P.y / Field->dx;
    coefs_x[2] = -coefs_x[3];
    coefs_x[1] = (1.0 - P.y) / Field->dx;
    coefs_x[0] = -coefs_x[1];

    double coefs_y[4];
    coefs_y[3] = P.x / Field->dy;
    coefs_y[2] = (1.0 - P.x) / Field->dy;
    coefs_y[1] = -coefs_y[3];
    coefs_y[0] = -coefs_y[2];

    Rank3 res;
    memset(&res, 0, sizeof(Rank3));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
        {
            ((Rank1 *)&res)[i].x += coefs_x[j] *((double *)elems[j])[i];
            ((Rank1 *)&res)[i].y += coefs_y[j] *((double *)elems[j])[i];
        }

    return res;
}

Rank3 christoffel_at_point(Rank2_field *Field, Rank1 P)
{
    Rank3 diff = diff_at_point(Field, P);

    Rank3 comb;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                ((double *)&comb)[4 * i + 2 * j + k] = 0.5 * (
                    ((double *)&diff)[4 * i + 2 * j + k] +
                    ((double *)&diff)[4 * i + 2 * k + j] -
                    ((double *)&diff)[4 * k + 2 * j + i]);

    Rank2 inv = value_at_point(Field, P);
    inv = rank2_inverse(&inv);

    Rank3 res;
    memset(&res, 0, sizeof(Rank3));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
            ((double *)&res)[4 * j + i] +=
                ((double *)&inv)[2 * j + k] * ((double *)&comb)[4 * k + i];

    return res;
}

void particle_init(Particle_2d *particle, const char* filename)
{
    FILE *const file = fopen(filename, "r");

    fscanf(file, "%lf %lf", &particle->pos.x, &particle->pos.y);
    fscanf(file, "%lf %lf", &particle->vel.x, &particle->vel.y);

    particle->time = 0;

    fclose(file);
}

void move_particle(Rank2_field *Field, Particle_2d *particle)
{
    Rank3 christ = christoffel_at_point(Field, particle->pos);

    const double ds = 1.0;
    const double *u = (const double *)&particle->vel;
    Rank1 du = {0, 0};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                ((double *)&du)[i] -=
                    ((double *)&christ)[4 * i + 2 * j + k] * u[k] * u[j];

    particle->pos.x += particle->vel.x * ds;
    if (particle->pos.x > 1.0)
    {
        particle->pos.x = 2.0 - particle->pos.x;
        particle->vel.x = - particle->vel.x;
    }
    if (particle->pos.x < 0.0)
    {
        particle->pos.x = - particle->pos.y;
        particle->vel.x = - particle->vel.x;
    }
    particle->pos.y += particle->vel.y * ds;
    if (particle->pos.y > 1.0)
        particle->pos.y -= 1.0;
    if (particle->pos.y < 0.0)
        particle->pos.y += 1.0;
    particle->vel.x += du.x * ds;
    particle->vel.y += du.y * ds;

    particle->time += ds;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        printf("Usage:\n");
        printf("\tgrav [input_file] [output_file]\n");
        return -1;
    }

    Rank2_field metric_field;
    Particle_2d tester;

    Rank2_field *Field = &metric_field;
    memset(Field, 0, sizeof(Rank2_field));
    rank2_field_alloc(Field, 10, 10);
    rank2_field_init(Field, calc_metric);

    Particle_2d *particle = &tester;
    memset(particle, 0, sizeof(Particle_2d));
    particle_init(particle, argv[1]);

    FILE *out = fopen(argv[2], "w");
    fprintf(out, "%lf\t%lf\t%lf\t%lf\t%lf\n", particle->time,
            particle->pos.x, particle->pos.y,
            particle->vel.x, particle->vel.y);
    for (int i = 0; i < 100; i++)
    {
        move_particle(Field, particle);
        fprintf(out, "%lf\t%lf\t%lf\t%lf\t%lf\n", particle->time,
                particle->pos.x, particle->pos.y,
                particle->vel.x, particle->vel.y);
    }
    fclose(out);

    delete_field(Field);

    return 0;
}
