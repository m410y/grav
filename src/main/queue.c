#include <queue.h>

#include <stdlib.h>

void queue_free(Queue *queue)
{
    Queue to_free;

    if (!(*queue))
        return;

    while ((*queue)->next != (*queue))
    {
        to_free = *queue;
        *queue = to_free->next;
        free(to_free);
    }

    free(*queue);
    *queue = NULL;
}

int queue_push(Queue *queue, void *data)
{
    Queue new_elem;

    new_elem = malloc(sizeof(Queue));
    if (!new_elem)
        return -1;

    new_elem->data = data;

    if (*queue)
    {
        new_elem->next = *queue;
        new_elem->prev = (*queue)->prev;
    }
    else
        new_elem->next = new_elem->prev = new_elem;

    *queue = new_elem;
    return 0;
}

void * queue_pop(Queue *queue)
{
    Queue to_free;
    void *data;

    if (*queue == NULL)
        return NULL;

    to_free = (*queue)->prev;
    data = to_free->data;

    if (to_free->prev != to_free)
    {
        (*queue)->prev = to_free->prev;
        to_free->prev->next = (*queue);
    }
    else
        *queue = NULL;

    free(to_free);
    return data;
}
