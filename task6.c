#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>

int compare(const void *a, const void *b)
{
    int int_a = *((int *)a);
    int int_b = *((int *)b);
    return (int_a > int_b) - (int_a < int_b);
}

void parallel_merge(int *array, int left, int mid, int right, int *temp)
{
    int length = right - left;
    if (length <= 10000)
    {
        int i = left, j = mid, k = left;
        while (i < mid && j < right)
        {
            if (array[i] <= array[j])
            {
                temp[k++] = array[i++];
            }
            else
            {
                temp[k++] = array[j++];
            }
        }
        while (i < mid)
        {
            temp[k++] = array[i++];
        }
        while (j < right)
        {
            temp[k++] = array[j++];
        }
        for (i = left; i < right; i++)
        {
            array[i] = temp[i];
        }
    }
    else
    {
        int mid1 = left + (mid - left) / 2;
        int mid2 = mid + (right - mid) / 2;

#pragma omp task shared(array, temp) firstprivate(left, mid1, mid)
        parallel_merge(array, left, mid1, mid, temp);

#pragma omp task shared(array, temp) firstprivate(mid, mid2, right)
        parallel_merge(array, mid, mid2, right, temp);

#pragma omp taskwait

        int i = left, j = mid, k = left;
        while (i < mid && j < right)
        {
            if (array[i] <= array[j])
            {
                temp[k++] = array[i++];
            }
            else
            {
                temp[k++] = array[j++];
            }
        }
        while (i < mid)
        {
            temp[k++] = array[i++];
        }
        while (j < right)
        {
            temp[k++] = array[j++];
        }
        for (i = left; i < right; i++)
        {
            array[i] = temp[i];
        }
    }
}

void parallel_merge_sort(int *array, int left, int right, int *temp, int depth)
{
    if (left >= right - 1)
    {
        return;
    }

    if (depth <= 0 || right - left <= 10000)
    {
        qsort(array + left, right - left, sizeof(int), compare);
    }
    else
    {
        int mid = left + (right - left) / 2;

#pragma omp task shared(array, temp) firstprivate(left, mid, depth)
        parallel_merge_sort(array, left, mid, temp, depth - 1);

#pragma omp task shared(array, temp) firstprivate(mid, right, depth)
        parallel_merge_sort(array, mid, right, temp, depth - 1);

#pragma omp taskwait

        parallel_merge(array, left, mid, right, temp);
    }
}

int main(int argc, char *argv[])
{
    int N = atoi(argv[1]);
    int p = atoi(argv[2]);

    int *A = (int *)malloc(N * sizeof(int));
    int *B = (int *)malloc(N * sizeof(int));
    int *temp = (int *)malloc(N * sizeof(int));

    if (A == NULL || B == NULL || temp == NULL)
    {
        free(A);
        free(B);
        free(temp);
        return 1;
    }

    srand(time(NULL));
    for (int i = 0; i < N; i++)
    {
        A[i] = rand();
        B[i] = A[i];
    }

    double start_time, end_time, parallel_time, qsort_time;

    omp_set_num_threads(p);

    int max_depth = 0;
    int threads = p;
    while (threads >>= 1)
        ++max_depth;

    start_time = omp_get_wtime();

#pragma omp parallel
    {
#pragma omp single
        {
            parallel_merge_sort(A, 0, N, temp, max_depth);
        }
    }

    end_time = omp_get_wtime();
    parallel_time = end_time - start_time;

    start_time = omp_get_wtime();
    qsort(B, N, sizeof(int), compare);
    end_time = omp_get_wtime();
    qsort_time = end_time - start_time;

    int correct = 1;
    for (int i = 0; i < N; i++)
    {
        if (A[i] != B[i])
        {
            correct = 0;
            printf("Несоответствие на позиции %d: %d != %d\n", i, A[i], B[i]);
            break;
        }
    }

    if (correct)
    {
        printf("Массив успешно отсортирован.\n");
    }
    else
    {
        printf("Ошибка сортировки массива.\n");
    }

    printf("Параллельная сортировка: %f секунд\n", parallel_time);
    printf("Стандартный qsort: %f секунд\n", qsort_time);

    if (parallel_time <= 1.05 * qsort_time)
    {
        printf("Многопоточный вариант соответствует требованию по времени.\n");
    }
    else
    {
        printf("Многопоточный вариант НЕ соответствует требованию по времени.\n");
    }

    free(A);
    free(B);
    free(temp);

    return 0;
}
