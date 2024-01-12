#include <list.h>

#include <stdlib.h>

void list_free(List *list)
{
    while (*list)
    {
        List to_free = *list;
        *list = to_free->next;
        free(to_free);
    }
}

int list_push(List *list, void *data)
{
    List new_elem = malloc(sizeof(List));
    if (!new_elem)
        return -1;

    new_elem->data = data;

    new_elem->next = *list;
    *list = new_elem;

    return 0;
}

void * list_pop(List *list)
{
    List to_free;
    void *data;

    if (*list == NULL)
        return NULL;

    to_free = *list;
    data = to_free->data;

    *list = to_free->next;

    free(to_free);

    return data;
}
