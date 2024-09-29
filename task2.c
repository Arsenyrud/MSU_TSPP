#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

typedef struct
{
    int *queue;
    int capacity;
    int head;
    int tail;
    int count;
    pthread_mutex_t mutex;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
} MyConcurrentQueue;

int init(MyConcurrentQueue *q, int capacity)
{
    q->capacity = capacity;
    q->count = 0;
    q->head = 0;
    q->tail = 0;
    q->queue = (int *)malloc(sizeof(int) * capacity);
    pthread_mutex_init(&q->mutex, NULL);
    pthread_cond_init(&q->not_empty, NULL);
    pthread_cond_init(&q->not_full, NULL);
    return 0;
}

void destroy_queue(MyConcurrentQueue *q)
{
    free(q->queue);
    pthread_mutex_destroy(&q->mutex);
    pthread_cond_destroy(&q->not_empty);
    pthread_cond_destroy(&q->not_full);
}

void put(MyConcurrentQueue *q, int value)
{
    pthread_mutex_lock(&q->mutex);
    while (q->count == q->capacity)
    {
        pthread_cond_wait(&q->not_full, &q->mutex);
    }
    q->queue[q->tail] = value;
    q->tail = (q->tail + 1) % q->capacity;
    q->count += 1;

    pthread_cond_signal(&q->not_empty);
    pthread_mutex_unlock(&q->mutex);
}

int get(MyConcurrentQueue *q)
{
    pthread_mutex_lock(&q->mutex);
    while (q->count == 0)
    {
        pthread_cond_wait(&q->not_empty, &q->mutex);
    }
    int value = q->queue[q->head];
    q->head = (q->head + 1) % q->capacity;
    q->count -= 1;
    pthread_cond_signal(&q->not_full);
    pthread_mutex_unlock(&q->mutex);
    return value;
}

typedef struct
{
    MyConcurrentQueue *q;
    int writer_id;
    int items_to_write;
} writer_args_t;

typedef struct
{
    MyConcurrentQueue *q;
    int reader_id;
    int items_to_read;
} reader_args_t;

void *writer(void *args)
{
    writer_args_t *wargs = (writer_args_t *)args;
    for (int i = 0; i < wargs->items_to_write; ++i)
    {
        int value = wargs->writer_id * 1000 + i;
        put(wargs->q, value);
        printf("Writer %d: put %d\n", wargs->writer_id, value);
        usleep(rand() % 100000);
    }
    return NULL;
}

void *reader(void *args)
{
    reader_args_t *rargs = (reader_args_t *)args;
    for (int i = 0; i < rargs->items_to_read; ++i)
    {
        int value = get(rargs->q);
        printf("Reader %d: got %d\n", rargs->reader_id, value);
        usleep(rand() % 150000);
    }
    return NULL;
}

int main()
{
    srand(time(NULL));

    int M = 2; // writers
    int N = 3; // readers
    int capacity = 5;
    int items_per_writer = 2;
    int total_items = M * items_per_writer;
    int items_per_reader = (total_items + N - 1) / N;

    pthread_t writers[M];
    pthread_t readers[N];

    MyConcurrentQueue q;
    init(&q, capacity);

    writer_args_t wargs[M];
    for (int i = 0; i < M; ++i)
    {
        wargs[i].q = &q;
        wargs[i].writer_id = i + 1;
        wargs[i].items_to_write = items_per_writer;
        pthread_create(&writers[i], NULL, writer, &wargs[i]);
    }

    reader_args_t rargs[N];
    for (int i = 0; i < N; ++i)
    {
        rargs[i].q = &q;
        rargs[i].reader_id = i + 1;
        if (i == N - 1)
        {
            rargs[i].items_to_read = total_items - items_per_writer * (M - 1);
        }
        else
        {
            rargs[i].items_to_read = items_per_reader;
        }
        pthread_create(&readers[i], NULL, reader, &rargs[i]);
    }

    for (int i = 0; i < M; ++i)
    {
        pthread_join(writers[i], NULL);
    }

    for (int i = 0; i < N; ++i)
    {
        pthread_join(readers[i], NULL);
    }

    destroy_queue(&q);

    return 0;
}