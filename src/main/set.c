#include <set.h>
#include <list.h>

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define START_SIZE 256
#define MAX_LOAD 1.0

/*
 * FNV-1
 */
static size_t hash(const void *data)
{
/*
 * 32 bit
 */
#if UINTPTR_MAX == 0xFFFFFFFF
    const int32_t = 0x811c9dc5 ;
    const int32_t prime = 0x01000193 ;
    const int bytes = sizeof(int32_t)/sizeof(unsigned char);
    int32_t res;
/*
 * 64 bit
 */
#else
    const int64_t bias = 0xcbf29ce484222325;
    const int64_t prime = 0x100000001b3;
    const int bytes = sizeof(int64_t)/sizeof(unsigned char);
    int64_t res;
#endif

    int i;
    unsigned char byte;

    res = bias;
    for (i = 0; i < bytes; i++)
    {
        byte = ((unsigned char *)&data)[i];
        res *= prime;
        res ^= byte;
    }

    return res;
}

static int set_init_size(Set *set, size_t size)
{
    set->size = size;
    set->elems = 0;
    set->table = malloc(set->size * sizeof(List *));
    if (!set->table)
        return -1;

    memset(set->table, 0, set->size * sizeof(List *));

    return 0;
}


static int set_add_nocheck(Set *set, void *data)
{
    list_push((List *)&set->table[hash(data) % set->size], data);
    set->elems++;

    return 0;
}

static int set_expand(Set *set)
{
    size_t i;
    List elem;
    struct set new_set;

    if (set_init_size(&new_set, set->size * 2) < 0)
        return -1;

    for (i = 0; i < set->size; i++)
        for (elem = set->table[i]; elem; elem = elem->next)
            set_add_nocheck(&new_set, elem->data);

    set_free(set);
    *set = new_set;

    return 0;
}

int set_init(Set *set)
{
    return set_init_size(set, START_SIZE);
}

void set_free(Set *set)
{
    size_t i;

    for (i = 0; i < set->size; i++)
        list_free((List *)&set->table[i]);

    free(set->table);
}

int set_contain(const Set *set, void *data)
{
    List elem;

    for (elem = set->table[hash(data) % set->size]; elem; elem = elem->next)
        if (elem->data == data)
            return 1;

    return 0;
}

int set_add(Set *set, void *data)
{
    if (set_contain(set, data))
        return 1;

    if ((double)set->elems / set->size > MAX_LOAD)
        if (set_expand(set) < 0)
            return -1;

    return set_add_nocheck(set, data);
}

int set_remove(Set *set, void *data)
{
    List list;

    list = set->table[hash(data) % set->size];
    if (list == NULL)
        return 1;

    if (list->data == data)
    {
        list_pop((List *)&set->table[hash(data) % set->size]);
        return 1;
    }

    for (; list->next; list = list->next)
        if (list->next->data == data)
        {
            list_pop(&list->next);
            return 1;
        }

    return 0;
}
