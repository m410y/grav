#include "tensor_deifnes.h"
#include <metric_field.h>
#include <set.h>
#include <queue.h>
#include <list.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

typedef struct face Face;
struct face {
    Face *near[3];
    Rank1 *vert[3];
    Rank2 appr[6];
};

static void gauss_sol_nocheck(int n, int m, double M[n][m])
{
    int i, j, k;
    double temp;

    for (i = 0; i < n; i++)
    {
        if (M[i][i] == 0.0)
        {
            for (k = i + 1; M[k][i] == 0.0; k++)
                ;

            for (int j = i; j < m; j++)
            {
                temp = M[i][j];
                M[i][j] = M[k][j];
                M[k][j] = temp;
            }
        }

        temp = 1.0 / M[i][i];
        for (j = i; j < m; j++)
            M[i][j] *= temp;

        for (k = i + 1; k < n; k++)
        {
            temp = M[k][i];
            for (j = i; j < m; j++)
                M[k][j] -= M[i][j] * temp;
        }
    }

    for (i = n - 1; i > 0; i--)
    {
        for (k = i - 1; k >= 0; k--)
        {
            temp = M[k][i];
            for (j = i; j < m; j++)
                M[k][j] -= M[i][j] * temp;
        }
    }
}

static void face_field_approx(Face *face, metric_func func)
{
    int i, j, k;
    double mat[6][10];
    Rank1 p[6];
    Rank2 g;

    for (i = 0; i < 3; i++)
    {
        j = (i + 1) % 3;
        k = (j + 1) % 3;

        p[i] = *(face->vert[i]);
        p[3 + k].x = 0.5 * (face->vert[j]->x + face->vert[k]->x);
        p[3 + k].y = 0.5 * (face->vert[j]->y + face->vert[k]->y);
    }

    for (i = 0; i < 6; i++)
    {
        g = func(p[i]);

        mat[i][0] = 1.0;
        mat[i][1] = p[i].x;
        mat[i][2] = p[i].y;
        mat[i][3] = p[i].x * p[i].x;
        mat[i][4] = p[i].x * p[i].y;
        mat[i][5] = p[i].y * p[i].y;

        for (j = 0; j < 4; j++)
            PTR_CAST(mat[i][6])[j] = PTR_CAST(g)[j];
    }

    gauss_sol_nocheck(6, 10, mat);

    for (i = 0; i < 6; i++)
        for (j = 0; j < 4; j++)
            PTR_CAST(face->appr[i])[j] = mat[i][6 + j];
}

static double triangle_area(Rank1 *p[3])
{
    double det = 0;
    int i, j, k;

    for (i = 0; i < 3; i++)
    {
        j = (i + 1) % 3;
        k = (j + 1) % 3;

        det += p[j]->x * p[k]->y;
        det -= p[k]->x * p[j]->y;
    }

    return det;
}

static int face_contain_point(Face *face, Rank1 *p)
{
    int i;
    Rank1 *temp;
    Rank1 *v[3];

    for (i = 0; i < 3; i++)
        v[i] = face->vert[i];

    if (triangle_area(v) < 0)
        for (i = 0; i < 3; i++)
            v[i] = face->vert[3 - i];

    for (i = 0; i < 3; i++)
    {
        temp = v[i];
        v[i] = p;
        if (triangle_area(v) < 0)
            return 0;

        v[i] = temp;
    }

    return 1;
}

static Face * field_face_at_point(Metric_field *field, Rank1 *P)
{
    int i;
    Set used;
    Queue queue;
    Face *face;

    set_init(&used);
    queue = NULL;

    set_add(&used, field->f);
    queue_push(&queue, field->f);
    while (queue)
    {
        face = queue_pop(&queue);
        for (i = 0; i < 3; i++)
        {
            if (!face->near[i])
                continue;

            if (set_add(&used, face->near[i]) == 0)
                queue_push(&queue, face->near[i]);

            if (face_contain_point(face, P))
                return field->f = face;
        }
    }

    set_free(&used);

    return NULL;
}

static Rank2 metric_at_point_nocheck(const Face *face, const Rank1 *P)
{
    int i;
    Rank2 res = {};

    for (i = 0; i < 4; i++)
    {
        PTR_CAST(res)[i] += PTR_CAST(face->appr[0])[i];
        PTR_CAST(res)[i] += PTR_CAST(face->appr[1])[i] * P->x;
        PTR_CAST(res)[i] += PTR_CAST(face->appr[2])[i] * P->y;
        PTR_CAST(res)[i] += PTR_CAST(face->appr[3])[i] * P->x * P->x;
        PTR_CAST(res)[i] += PTR_CAST(face->appr[4])[i] * P->x * P->y;
        PTR_CAST(res)[i] += PTR_CAST(face->appr[5])[i] * P->y * P->y;
    }

    return res;
}

static Rank3 diff_at_point_nocheck(const Face *face, const Rank1 *P)
{
    int i;
    Rank3 res = {};

    for (i = 0; i < 4; i++)
    {
        PTR_CAST(res)[i] += PTR_CAST(face->appr[1])[i];
        PTR_CAST(res)[i] += PTR_CAST(face->appr[3])[i] * 2.0 * P->x;
        PTR_CAST(res)[i] += PTR_CAST(face->appr[4])[i] * P->y;

        PTR_CAST(res)[4 + i] += PTR_CAST(face->appr[2])[i];
        PTR_CAST(res)[4 + i] += PTR_CAST(face->appr[4])[i] * P->x;
        PTR_CAST(res)[4 + i] += PTR_CAST(face->appr[5])[i] * 2.0 * P->y;
    }

    return res;
}

static Rank4 second_at_point_nocheck(const Face *face)
{
    int i;
    Rank4 res = {};

    for (i = 0; i < 4; i++)
    {
        PTR_CAST(res)[4*0 + i] += PTR_CAST(face->appr[3])[i] * 2.0;
        PTR_CAST(res)[4*1 + i] += PTR_CAST(face->appr[4])[i];
        PTR_CAST(res)[4*2 + i] += PTR_CAST(face->appr[4])[i];
        PTR_CAST(res)[4*3 + i] += PTR_CAST(face->appr[5])[i] * 2.0;
    }

    return res;
}

static Rank2 inverse(const Rank2 *M)
{
    const double det = M->x.x * M->y.y - M->x.y * M->y.x;
    const Rank2 res = {
        .x.x = M->y.y / det,
        .x.y = -M->x.y / det,
        .y.x = -M->y.x / det,
        .y.y = M->x.x / det,
    };

    return res;
}

int field_alloc(Metric_field *field, int Nx, int Ny)
{
    Rank1 *verts;
    Face **faces;
    int i, j, k;

    verts = malloc((Nx + 1) * (Ny + 1) * sizeof(Rank1));
    if (!verts)
        return -1;

    for (i = 0; i <= Ny; i++)
        for (j = 0; j <= Nx; j++)
        {
            k = (Nx + 1) * i + j;
            verts[k].x = (double)j / Nx;
            verts[k].y = (double)i / Ny;
        }

    faces = malloc(2 * Nx * Ny * sizeof(Face));
    if (!faces)
    {
        free(verts);
        return -1;
    }

    for (i = 0; i < Ny; i++)
        for (j = 0; j < Nx; j++)
        {
            k = 2 * (i * Nx + j);

            faces[k] = malloc(sizeof(Face));
            faces[k + 1] = malloc(sizeof(Face));
            if (!faces[k] || !faces[k + 1])
            {
                i = Ny;
                j = faces[k] ? k + 1 : k;
                break;
            }

            memset(faces[k], 0, sizeof(Face));
            memset(faces[k + 1], 0, sizeof(Face));

            faces[k]->near[0] = faces[k + 1];
            faces[k + 1]->near[0] = faces[k];

            if (j != 0)
            {
                faces[k + 1]->near[1] = faces[k - 2];
                faces[k - 2]->near[1] = faces[k + 1];
            }

            if (i != 0)
            {
                faces[k]->near[2] = faces[k - 2 * Nx + 1];
                faces[k - 2 * Nx + 1]->near[2] = faces[k];
            }

            faces[k]->vert[0] = &verts[(Nx + 1) * i + j];
            faces[k]->vert[1] = &verts[(Nx + 1) * i + (j + 1)];
            faces[k]->vert[2] = &verts[(Nx + 1) * (i + 1) + (j + 1)];

            faces[k + 1]->vert[0] = &verts[(Nx + 1) * i + j];
            faces[k + 1]->vert[1] = &verts[(Nx + 1) * (i + 1) + (j + 1)];
            faces[k + 1]->vert[2] = &verts[(Nx + 1) * (i + 1) + j];
        }

    if (!faces[k] || !faces[k + 1])
    {
        for (i = j; i >= 0; i--)
            free(faces[i]);

        free(faces);
        free(verts);
        return -1;
    }

    field->f = faces[Nx * Ny];
    field->vertices = verts;

    free(faces);
    return 0;
}

void field_delete(Metric_field *field)
{
    int i;
    Set used;
    List list;
    Face *face;

    set_init(&used);
    list = NULL;

    set_add(&used, field->f);
    list_push(&list, field->f);
    while (list)
    {
        face = list_pop(&list);
        for (i = 0; i < 3; i++)
        {
            if (!face->near[i])
                continue;

            if (set_add(&used, face->near[i]) == 0)
                list_push(&list, face->near[i]);
        }

        free(face);
    }

    set_free(&used);
    list_free(&list);

    free(field->vertices);

    field->f = NULL;
    field->vertices = NULL;
}

void field_func_init(Metric_field *field, metric_func func)
{
    int i;
    Set used;
    List list;
    Face *face;

    set_init(&used);
    list = NULL;

    set_add(&used, field->f);
    list_push(&list, field->f);
    while (list)
    {
        face = list_pop(&list);
        for (i = 0; i < 3; i++)
        {
            if (!face->near[i])
                continue;

            if (set_add(&used, face->near[i]) == 0)
                list_push(&list, face->near[i]);

            face_field_approx(face, func);
        }
    }

    set_free(&used);
    list_free(&list);
}

Rank2 field_metric_at_point(Metric_field *field, Rank1 P)
{
    const Face *face = field_face_at_point(field, &P);
    return metric_at_point_nocheck(face, &P);
}

Rank3 field_christoffel_at_point(Metric_field *field, Rank1 P)
{
    Face *face;
    int i, j, k;
    Rank2 inv;
    Rank3 comb, diff;
    Rank3 res = {};

    face = field_face_at_point(field, &P);

    diff = diff_at_point_nocheck(face, &P);
    inv = metric_at_point_nocheck(face, &P);
    inv = inverse(&inv);

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                PTR_CAST(comb)[4*i + 2*j + k] = 0.5 * (
                    PTR_CAST(diff)[4*k + 2*i + j] +
                    PTR_CAST(diff)[4*j + 2*i + k] -
                    PTR_CAST(diff)[4*i + 2*j + k]);

    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                PTR_CAST(res)[4 * j + i] +=
                    PTR_CAST(inv)[2*j + k] * PTR_CAST(comb)[4*k + i];

    return res;
}

Rank4 field_riemann_at_point(Metric_field *field, Rank1 P)
{
    int i, j, k, l, m, idx;
    Face *face;
    Rank2 inv;
    Rank3 christ;
    Rank4 diff;
    double temp;
    Rank4 res = {};

    face = field_face_at_point(field, &P);

    christ = field_christoffel_at_point(field, P);
    diff = second_at_point_nocheck(face);
    inv = metric_at_point_nocheck(face, &P);
    inv = inverse(&inv);

    temp = diff.x.y.x.y - 0.5*(diff.x.x.y.y + diff.y.y.x.x);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
            {
                m = (j + 1) % 2;
                l = (k + 1) % 2;

                idx = 8*i + 4*j + 2*k + l;
                PTR_CAST(res)[idx] =
                    (j ? temp : -temp) * PTR_CAST(inv)[2*i + m];

                for (m = 0; m < 2; m++)
                {
                    PTR_CAST(res)[idx] +=
                        PTR_CAST(christ)[4*i + 2*k + m] *
                        PTR_CAST(christ)[4*m + 2*j + l] -
                        PTR_CAST(christ)[4*i + 2*l + m] *
                        PTR_CAST(christ)[4*m + 2*j + k];
                }
            }

    return res;
}

Rank2 field_ricci_at_point(Metric_field *field, Rank1 P)
{
    int i, j, k;
    Rank4 riemann;
    Rank2 res = {};

    riemann = field_riemann_at_point(field, P);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                PTR_CAST(res)[2*i + j] += PTR_CAST(riemann)[8*k + 4*i + 2*k + j];

    return res;
}

Rank0 field_curavture_at_point(Metric_field *field, Rank1 P)
{
    int i;
    Rank2 ricci, inv;
    Rank0 res = 0;

    ricci = field_ricci_at_point(field, P);
    inv = field_metric_at_point(field, P);
    inv = inverse(&inv);

    for (i = 0; i < 4; i++)
        res += PTR_CAST(ricci)[i] * PTR_CAST(inv)[i];

    return res;
}
