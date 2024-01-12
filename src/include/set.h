#ifndef SET_H
#define SET_H

#include <stddef.h>

typedef struct set Set;
struct set {
    void **table;
    size_t size;
    size_t elems;
};

int set_init(Set *set);
void set_free(Set *set);
int set_contain(const Set *set, void *elem);
int set_add(Set *set, void *elem);
int set_remove(Set *set, void *elem);

#endif
