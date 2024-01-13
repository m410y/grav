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
    Rank2 metric;
    Rank3 d_metric;
};

static void face_field_approx(Face *face, metric_func func)
{
    int i, j, k;
    double areas[3];
    double y_diff[3];
    double x_diff[3];
    double det;
    Rank2 vals[3];
    const Rank1 *p[3] = {
        face->vert[0],
        face->vert[1],
        face->vert[2],
    };

    det = 0;
    for (i = 0; i < 3; i++)
    {
        j = (i + 1) % 3;
        k = (j + 1) % 3;

        areas[i] = p[j]->x * p[k]->y - p[k]->x * p[j]->y;
        y_diff[i] = p[j]->y - p[k]->y;
        x_diff[i] = p[k]->x - p[j]->x;

        det += areas[i];

        vals[i] = func(*p[i]);
    }

    memset(&face->metric, 0, sizeof(face->metric));
    memset(&face->d_metric, 0, sizeof(face->d_metric));

    for (i = 0; i < 3; i++)
    {
        areas[i] /= det;
        y_diff[i] /= det;
        x_diff[i] /= det;

        for (j = 0; j < 4; j++)
        {
            PTR_CAST(face->metric)[j] +=
                areas[i] * PTR_CAST(vals[i])[j];

            PTR_CAST(face->d_metric.x)[j] +=
                y_diff[i] * PTR_CAST(vals[i])[j];

            PTR_CAST(face->d_metric.y)[j] +=
                x_diff[i] * PTR_CAST(vals[i])[j];
        }
    }
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
    Rank2 res;

    res = face->metric;
    for (i = 0; i < 4; i++)
    {
        PTR_CAST(res)[i] += PTR_CAST(face->d_metric.x)[i] * P->x;
        PTR_CAST(res)[i] += PTR_CAST(face->d_metric.y)[i] * P->y;
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
    Rank3 comb, res;

    face = field_face_at_point(field, &P);
    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
                PTR_CAST(comb)[4 * i + 2 * j + k] = 0.5 * (
                    PTR_CAST(face->d_metric)[4 * k + 2 * i + j] +
                    PTR_CAST(face->d_metric)[4 * j + 2 * i + k] -
                    PTR_CAST(face->d_metric)[4 * i + 2 * j + k]);

    inv = metric_at_point_nocheck(face, &P);
    inv = inverse(&inv);
    memset(&res, 0, sizeof(Rank3));
    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            for (k = 0; k < 2; k++)
            PTR_CAST(res)[4 * j + i] +=
                PTR_CAST(inv)[2 * j + k] * PTR_CAST(comb)[4 * k + i];

    return res;
}
