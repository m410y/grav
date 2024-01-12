#ifndef QUEUE_H
#define QUEUE_H

typedef struct queue * Queue;
struct queue {
    Queue next;
    Queue prev;
    void *data;
};

void queue_free(Queue *queue);
int queue_push(Queue *queue, void *data);
void * queue_pop(Queue *queue);

#endif
