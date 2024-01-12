#include <particle.h>

static Rank1 particle_acceleration(Metric_field *field, const Rank1 *P, const Rank1 *V)
{
    int i, j, k;
    Rank3 christ;
    Rank1 du;
    
    christ = christoffel_at_point(field, *P);
    du.x = 0.0;
    du.y = 0.0;
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                PTR_CAST(du)[i] -=
                    PTR_CAST(christ)[4 * i + 2 * j + k] *
                    PTR_CAST(V)[j] * PTR_CAST(V)[k];

    return du;
}

void move_particle(Particle *particle, Metric_field *field, double ds)
{
    const double coeff1[4] = {
        0.0, 0.5, 0.5, 1.0,
    };
    const double coeff2[4] = {
        1.0 / 6, 2.0 / 6, 2.0 / 6, 1.0 / 6,
    };
    int i, j;
    Rank1 P, V;
    Rank1 dP[4], dV[4];

    for (i = 0; i < 4; i++)
    {
        P = particle->pos;
        V = particle->vel;
        if (coeff1[i] != 0.0)
            for (j = 0; j < 2; j++)
            {
                PTR_CAST(P)[j] += PTR_CAST(dP[i])[j] * coeff1[i] * ds;
                PTR_CAST(V)[j] += PTR_CAST(dV[i])[j] * coeff1[i] * ds;
            }

        dP[i] = V;
        dV[i] = particle_acceleration(field, &P, &V);
    }

    P = particle->pos;
    V = particle->vel;
    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
        {
            PTR_CAST(P)[j] += PTR_CAST(dP[i])[j] * coeff2[i] * ds;
            PTR_CAST(V)[j] += PTR_CAST(dV[i])[j] * coeff2[i] * ds;
        }

    particle->pos = P;
    particle->vel = V;
    particle->time += ds;
}
