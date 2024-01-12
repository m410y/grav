#ifdef UNIT_TESTS

#include <list.h>

#include <check.h>
#include <stdlib.h>


START_TEST (heap_test)
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

START_TEST (stack_test)
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

START_TEST (pop_null_test)
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

START_TEST (free_test)
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

Suite * list_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("List");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, heap_test);
    tcase_add_test(tc_core, stack_test);
    tcase_add_test(tc_core, pop_null_test);
    tcase_add_test(tc_core, free_test);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int number_failed;
    Suite *s;
    SRunner *sr;

    s = list_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    number_failed = srunner_ntests_failed(sr);

    srunner_free(sr);
    return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#endif
