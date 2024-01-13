#include "metric_field.h"
#ifdef UNIT_TESTS

#include <check.h>
#include <stdlib.h>

#include <list.h>
#include <queue.h>
#include <set.h>
#include <metric_field.h>
#include <particle.h>

START_TEST (list_heap_test)
{
    List *plist;

    plist = malloc(sizeof(List));
    ck_assert_ptr_nonnull(plist);
    *plist = NULL;
    ck_assert_ptr_null(*plist);
    list_push(plist, (void *)1);
    ck_assert_ptr_nonnull(*plist);
    ck_assert_ptr_eq(list_pop(plist), (void *)1);
    ck_assert_ptr_null(*plist);
    free(plist);
}
END_TEST

START_TEST (list_stack_test)
{
    List list;

    list = NULL;
    ck_assert_ptr_nonnull(&list);
    ck_assert_ptr_null(list);
    list_push(&list, (void *)1);
    ck_assert_ptr_nonnull(list);
    ck_assert_ptr_eq(list_pop(&list), (void *)1);
    ck_assert_ptr_null(list);
}
END_TEST

START_TEST (list_pop_null_test)
{
    List list;

    list = NULL;
    list_push(&list, (void *)1);
    list_push(&list, (void *)2);
    list_push(&list, (void *)3);

    ck_assert_ptr_eq(list_pop(&list), (void *)3);
    ck_assert_ptr_eq(list_pop(&list), (void *)2);
    ck_assert_ptr_eq(list_pop(&list), (void *)1);

    ck_assert_ptr_eq(list_pop(&list), NULL);
    ck_assert_ptr_eq(list_pop(&list), NULL);

    ck_assert_ptr_null(list);
}
END_TEST

START_TEST (list_free_test)
{
    List list;

    list = NULL;
    list_push(&list, (void *)1);
    list_push(&list, (void *)2);
    list_push(&list, (void *)3);

    list_free(&list);
    ck_assert_ptr_null(list);
}
END_TEST

START_TEST (queue_heap_test)
{
    Queue *pqueue;

    pqueue = malloc(sizeof(Queue));
    ck_assert_ptr_nonnull(pqueue);
    *pqueue = NULL;
    ck_assert_ptr_null(*pqueue);
    queue_push(pqueue, (void *)1);
    ck_assert_ptr_nonnull(*pqueue);
    ck_assert_ptr_eq(queue_pop(pqueue), (void *)1);
    ck_assert_ptr_null(*pqueue);
    free(pqueue);
}
END_TEST

START_TEST (queue_stack_test)
{
    Queue queue;

    queue = NULL;
    ck_assert_ptr_nonnull(&queue);
    ck_assert_ptr_null(queue);
    queue_push(&queue, (void *)1);
    ck_assert_ptr_nonnull(queue);
    ck_assert_ptr_eq(queue_pop(&queue), (void *)1);
    ck_assert_ptr_null(queue);
}
END_TEST

START_TEST (queue_pop_null_test)
{
    Queue queue;

    queue = NULL;
    queue_push(&queue, (void *)1);
    queue_push(&queue, (void *)2);
    queue_push(&queue, (void *)3);

    ck_assert_ptr_eq(queue_pop(&queue), (void *)1);
    ck_assert_ptr_eq(queue_pop(&queue), (void *)2);
    ck_assert_ptr_eq(queue_pop(&queue), (void *)3);

    ck_assert_ptr_eq(queue_pop(&queue), NULL);
    ck_assert_ptr_eq(queue_pop(&queue), NULL);

    ck_assert_ptr_null(queue);
}
END_TEST

START_TEST (queue_free_test)
{
    Queue queue;

    queue = NULL;
    queue_push(&queue, (void *)1);
    queue_push(&queue, (void *)2);
    queue_push(&queue, (void *)3);

    queue_free(&queue);
    ck_assert_ptr_null(queue);
}
END_TEST

START_TEST (set_heap_test)
{
    Set *pset;

    pset = malloc(sizeof(Set));
    ck_assert_ptr_nonnull(pset);
    set_init(pset);
    ck_assert_ptr_nonnull(pset->table);
    ck_assert_uint_ne(pset->size, 0);
    ck_assert_uint_eq(pset->elems, 0);
    set_free(pset);
    free(pset);
}
END_TEST

START_TEST (set_stack_test)
{
    Set set;

    set_init(&set);
    ck_assert_ptr_nonnull(set.table);
    ck_assert_uint_ne(set.size, 0);
    ck_assert_uint_eq(set.elems, 0);
    set_free(&set);
}
END_TEST

START_TEST (set_add_remove_test)
{
    Set set;

    set_init(&set);
    set_add(&set, (void *)1);
    set_remove(&set, (void *)1);
    set_free(&set);
}
END_TEST

START_TEST (set_contain_test)
{
    Set set;

    set_init(&set);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    ck_assert_int_eq(set_contain(&set, (void *)2), 0);
    set_add(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 1);
    ck_assert_int_eq(set_contain(&set, (void *)2), 0);
    set_add(&set, (void *)2);
    ck_assert_int_eq(set_contain(&set, (void *)1), 1);
    ck_assert_int_eq(set_contain(&set, (void *)2), 1);
    set_remove(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    ck_assert_int_eq(set_contain(&set, (void *)2), 1);
    set_remove(&set, (void *)2);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    ck_assert_int_eq(set_contain(&set, (void *)2), 0);
    set_free(&set);
}
END_TEST

START_TEST (set_add_more_test)
{
    Set set;

    set_init(&set);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    set_add(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 1);
    set_add(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 1);
    set_remove(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    set_free(&set);
}
END_TEST

START_TEST (set_remove_more_test)
{
    Set set;

    set_init(&set);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    set_add(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 1);
    set_remove(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    set_remove(&set, (void *)1);
    ck_assert_int_eq(set_contain(&set, (void *)1), 0);
    set_free(&set);
}
END_TEST

START_TEST (set_free_nonempty_test)
{
    Set set;

    set_init(&set);
    set_add(&set, (void *)1);
    set_add(&set, (void *)1);
    set_add(&set, (void *)2);
    set_add(&set, (void *)3);
    set_free(&set);
}
END_TEST

START_TEST (set_remove_empty_test)
{
    Set set;

    set_init(&set);
    set_remove(&set, (void *)1);
    set_remove(&set, (void *)2);
    set_remove(&set, (void *)3);
    set_add(&set, (void *)1);
    set_free(&set);
}
END_TEST

static Rank2 metric(Rank1 p)
{
    const Rank2 res = {
        {1.0 + p.x * p.x, p.x * p.y},
        {p.x * p.y, 1.0 + p.y * p.y},
    };
    return res;
}

START_TEST (mfield_heap_test)
{
    Metric_field *pfield;

    pfield = malloc(sizeof(Metric_field));
    ck_assert_ptr_nonnull(pfield);
    field_alloc(pfield, 10, 10);
    ck_assert_ptr_nonnull(pfield->f);
    ck_assert_ptr_nonnull(pfield->vertices);
    field_delete(pfield);
    free(pfield);
}
END_TEST

START_TEST (mfield_stack_test)
{
    Metric_field field;

    field_alloc(&field, 10, 10);
    ck_assert_ptr_nonnull(field.f);
    ck_assert_ptr_nonnull(field.vertices);
    field_delete(&field);
}
END_TEST

START_TEST (mfield_init_test)
{
    Metric_field field;

    field_alloc(&field, 10, 10);
    field_func_init(&field, &metric);
    field_delete(&field);
}
END_TEST

START_TEST (mfield_metric_test)
{
    const double eps = 1e-6;
    const Rank1 p = {
        0.5, 0.5,
    };
    const Rank2 ans = {
        {1.0 + p.x*p.x, p.x*p.y},
        {p.x*p.y, 1.0 + p.y*p.y},
    };
    int i;
    Metric_field field;
    Rank2 res;

    field_alloc(&field, 10, 10);
    field_func_init(&field, &metric);
    res = field_metric_at_point(&field, p);
    for (i = 0; i < 4; i++)
    {
        ck_assert_double_ge(PTR_CAST(res)[i], PTR_CAST(ans)[i] - eps);
        ck_assert_double_le(PTR_CAST(res)[i], PTR_CAST(ans)[i] + eps);
    }
    field_delete(&field);
}
END_TEST

START_TEST (mfield_christ_test)
{
    const double eps = 1e-2;
    const Rank1 p = {
        0.5, 0.5,
    };
    const Rank3 ans = {
        {{p.x/(1.0 + p.x*p.x + p.y*p.y), 0.0},
         {0.0, p.x/(1.0 + p.x*p.x + p.y*p.y)}},

        {{p.y/(1.0 + p.x*p.x + p.y*p.y), 0.0},
         {0.0, p.y/(1.0 + p.x*p.x + p.y*p.y)}},
    };
    int i;
    Metric_field field;
    Rank3 res;

    field_alloc(&field, 101, 101);
    field_func_init(&field, &metric);
    res = field_christoffel_at_point(&field, p);
    for (i = 0; i < 8; i++)
    {
        ck_assert_double_ge(PTR_CAST(res)[i], PTR_CAST(ans)[i] - eps);
        ck_assert_double_le(PTR_CAST(res)[i], PTR_CAST(ans)[i] + eps);
    }
    field_delete(&field);
}
END_TEST

START_TEST (particle_move_test)
{
    const double eps = 1e-4;
    int i;
    Metric_field field;
    Particle particle = {
        .pos = {0.5, 0.5},
        .vel = {0.1, 0.1},
        .time = 0.0,
    };
    const Particle ans = {
        .pos = {0.5 + 0.1, 0.5 + 0.1},
        .vel = {0.1 - 0.02/3, 0.1 - 0.02/3},
        .time = 1.0,
    };

    field_alloc(&field, 101, 101);
    field_func_init(&field, &metric);
    move_particle(&particle, &field, 1.0);
    for (i = 0; i < 5; i++)
    {
        ck_assert_double_ge(PTR_CAST(particle)[i], PTR_CAST(ans)[i] - eps);
        ck_assert_double_le(PTR_CAST(particle)[i], PTR_CAST(ans)[i] + eps);
    }
    field_delete(&field);
}
END_TEST

Suite * grav_suite(void)
{
    Suite *s;
    TCase *tc_list, *tc_queue, *tc_set, *tc_mfield, *tc_particle;

    s = suite_create("Grav");

    tc_list = tcase_create("List");
    tcase_add_test(tc_list, list_heap_test);
    tcase_add_test(tc_list, list_stack_test);
    tcase_add_test(tc_list, list_pop_null_test);
    tcase_add_test(tc_list, list_free_test);

    tc_queue = tcase_create("Queue");
    tcase_add_test(tc_queue, queue_heap_test);
    tcase_add_test(tc_queue, queue_stack_test);
    tcase_add_test(tc_queue, queue_pop_null_test);
    tcase_add_test(tc_queue, queue_free_test);

    tc_set = tcase_create("Set");
    tcase_add_test(tc_set, set_heap_test);
    tcase_add_test(tc_set, set_stack_test);
    tcase_add_test(tc_set, set_add_remove_test);
    tcase_add_test(tc_set, set_contain_test);
    tcase_add_test(tc_set, set_add_more_test);
    tcase_add_test(tc_set, set_remove_more_test);
    tcase_add_test(tc_set, set_free_nonempty_test);
    tcase_add_test(tc_set, set_remove_empty_test);

    tc_mfield = tcase_create("Metric_field");
    tcase_add_test(tc_mfield, mfield_heap_test);
    tcase_add_test(tc_mfield, mfield_stack_test);
    tcase_add_test(tc_mfield, mfield_init_test);
    tcase_add_test(tc_mfield, mfield_metric_test);
    tcase_add_test(tc_mfield, mfield_christ_test);

    tc_particle = tcase_create("Metric_field");
    tcase_add_test(tc_particle, particle_move_test);

    suite_add_tcase(s, tc_list);
    suite_add_tcase(s, tc_queue);
    suite_add_tcase(s, tc_set);
    suite_add_tcase(s, tc_mfield);
    suite_add_tcase(s, tc_particle);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = grav_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);

    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#endif
