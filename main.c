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

typedef struct rank2_field Rank2_field;
struct rank2_field {
    Rank2 *F;
    size_t Nx, Ny;
    double dx, dy;
    Rank1 (*conn)(Rank2_field *, Rank1 *, Rank2 **);
};

typedef struct particle_2d Particle_2d;
struct particle_2d {
    Rank1 pos;
    Rank1 vel;
    double time;
};

void field_alloc(Rank2_field *Field, size_t Nx, size_t Ny)
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
    Field->conn = NULL;
}

void field_func_init(Rank2_field *Field, Rank2 (*func)(Rank1))
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

Rank2 inverse(Rank2 *M)
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

Rank1 torus_conn(Rank2_field *Field, Rank1 *P, Rank2 *elems[4])
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

    if (elems)
    {
        elems[0] = &Field->F[i * Field->Nx + j];
        elems[1] = &Field->F[i * Field->Nx + jp];
        elems[2] = &Field->F[ip * Field->Nx + j];
        elems[3] = &Field->F[ip * Field->Nx + jp];
    }

    return Pp;
}

Rank2 value_at_point(Rank2_field *Field, Rank1 P)
{
    Rank2 *elems[4];
    P = Field->conn(Field, &P, elems);

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
    P = Field->conn(Field, &P, elems);

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
    inv = inverse(&inv);

    Rank3 res;
    memset(&res, 0, sizeof(Rank3));
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
            ((double *)&res)[4 * j + i] +=
                ((double *)&inv)[2 * j + k] * ((double *)&comb)[4 * k + i];

    return res;
}

void particle_init(Particle_2d *particle, FILE *const file)
{
    fscanf(file, "%lf %lf", &particle->pos.x, &particle->pos.y);
    fscanf(file, "%lf %lf", &particle->vel.x, &particle->vel.y);

    particle->time = 0;
}

Rank1 particle_acceleration(Rank2_field *Field, const Rank1 *P, const Rank1 *V)
{
    Rank3 christ = christoffel_at_point(Field, *P);
    const double *u = (const double *)V;
    Rank1 du = {0, 0};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                ((double *)&du)[i] -=
                    ((double *)&christ)[4 * i + 2 * j + k] * u[j] * u[k];

    return du;
}

void move_particle(Rank2_field *Field, Particle_2d *particle, double ds)
{
    Rank1 P = particle->pos;
    Rank1 V = particle->vel;

    Rank1 dP_1 = V;
    Rank1 dV_1 = particle_acceleration(Field, &P, &V);

    P.x = particle->pos.x + dP_1.x * ds / 2;
    P.y = particle->pos.y + dP_1.y * ds / 2;
    V.x = particle->vel.x + dV_1.x * ds / 2;
    V.y = particle->vel.y + dV_1.y * ds / 2;

    Rank1 dP_2 = V;
    Rank1 dV_2 = particle_acceleration(Field, &P, &V);

    P.x = particle->pos.x + dP_2.x * ds / 2;
    P.y = particle->pos.y + dP_2.y * ds / 2;
    V.x = particle->vel.x + dV_2.x * ds / 2;
    V.y = particle->vel.y + dV_2.y * ds / 2;

    Rank1 dP_3 = V;
    Rank1 dV_3 = particle_acceleration(Field, &P, &V);

    P.x = particle->pos.x + dP_3.x * ds;
    P.y = particle->pos.y + dP_3.y * ds;
    V.x = particle->vel.x + dV_3.x * ds;
    V.y = particle->vel.y + dV_3.y * ds;

    Rank1 dP_4 = V;
    Rank1 dV_4 = particle_acceleration(Field, &P, &V);

    P.x = particle->pos.x + (dP_1.x + 2 * dP_2.x + 2 * dP_3.x + dP_4.x) * ds / 6;
    P.y = particle->pos.y + (dP_1.y + 2 * dP_2.y + 2 * dP_3.y + dP_4.y) * ds / 6;
    V.x = particle->vel.x + (dV_1.x + 2 * dV_2.x + 2 * dV_3.x + dV_4.x) * ds / 6;
    V.y = particle->vel.y + (dV_1.y + 2 * dV_2.y + 2 * dV_3.y + dV_4.y) * ds / 6;

    P.x -= floor(P.x);
    P.y -= floor(P.y);

    particle->pos = P;
    particle->vel = V;

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
    int iterations;
    double step;
    int Nx, Ny;
    FILE *in = fopen(argv[1], "r");
    FILE *out = fopen(argv[2], "w");;

    fscanf(in, "%d", &iterations);
    fscanf(in, "%lf", &step);
    fscanf(in, "%d%d", &Nx, &Ny);

    Rank2_field *Field = &metric_field;
    memset(Field, 0, sizeof(Rank2_field));
    field_alloc(Field, Nx, Ny);
    field_func_init(Field, calc_metric);
    Field->conn = torus_conn;

    Particle_2d *particle = &tester;
    memset(particle, 0, sizeof(Particle_2d));
    particle_init(particle, in);

    fclose(in);

    fprintf(out, "%lf\t%lf\t%lf\t%lf\t%lf\n", particle->time,
            particle->pos.x, particle->pos.y,
            particle->vel.x, particle->vel.y);
    for (int i = 0; i < iterations; i++)
    {
        move_particle(Field, particle, step);
        fprintf(out, "%lf\t%lf\t%lf\t%lf\t%lf\n", particle->time,
                particle->pos.x, particle->pos.y,
                particle->vel.x, particle->vel.y);
    }

    fclose(out);

    delete_field(Field);

    return 0;
}
