#ifndef LIST_H
#define LIST_H

typedef struct list * List;
struct list {
    List next;
    void *data;
};

void list_free(List *list);
int list_push(List *list, void *data);
void * list_pop(List *list);

#endif
