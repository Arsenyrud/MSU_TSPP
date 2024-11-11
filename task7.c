#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 1024
#define MAX_ITER 1000
#define TOL 1e-6

int main(int argc, char *argv[]) {
    int rank, size;
    int i, j, iter;
    int n_local;
    double **f_old, **f_new;
    double diff, global_diff;
    double start_time, end_time;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (N % size != 0) {
        if (rank == 0) {
            fprintf(stderr, "N (%d) не делится на количество процессов (%d)\n", N, size);
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    n_local = N / size;
    f_old = (double **)malloc((n_local + 2) * sizeof(double *));
    f_new = (double **)malloc((n_local + 2) * sizeof(double *));
    for (i = 0; i < n_local + 2; i++) {
        f_old[i] = (double *)malloc(N * sizeof(double));
        f_new[i] = (double *)malloc(N * sizeof(double));
    }

    srand(rank + 1);
    for (i = 1; i <= n_local; i++) {
        for (j = 0; j < N; j++) {
            f_old[i][j] = ((double)rand()) / RAND_MAX;
            f_new[i][j] = f_old[i][j];
        }
    }

    if (rank == 0) {
        for (j = 0; j < N; j++) {
            f_old[0][j] = 0.0;
            f_new[0][j] = 0.0;
        }
    }
    if (rank == size - 1) {
        for (j = 0; j < N; j++) {
            f_old[n_local + 1][j] = 0.0;
            f_new[n_local + 1][j] = 0.0;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    start_time = MPI_Wtime();

    for (iter = 0; iter < MAX_ITER; iter++) {
       if (rank > 0) {
            MPI_Send(f_old[1], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(f_old[0], N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            for (j = 0; j < N; j++) {
                f_old[0][j] = 0.0;
            }
        }

        if (rank < size - 1) {
            MPI_Send(f_old[n_local], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(f_old[n_local + 1], N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } else {
            for (j = 0; j < N; j++) {
                f_old[n_local + 1][j] = 0.0;
            }
        }

        double local_diff = 0.0;
        for (i = 1; i <= n_local; i++) {
            for (j = 1; j < N - 1; j++) {
                f_new[i][j] = 0.25 * (f_old[i + 1][j] + f_old[i - 1][j] +
                                      f_old[i][j + 1] + f_old[i][j - 1]);
                double temp_diff = fabs(f_new[i][j] - f_old[i][j]);
                if (temp_diff > local_diff) {
                    local_diff = temp_diff;
                }
            }
        }

        MPI_Allreduce(&local_diff, &diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        for (i = 1; i <= n_local; i++) {
            for (j = 1; j < N - 1; j++) {
                f_old[i][j] = f_new[i][j];
            }
        }

        if (iter % 100 == 0 && rank == 0) {
            printf("Итерация %d, Макс разница: %f\n", iter, diff);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); 
    end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Максимальная разница после %d итераций: %f\n", MAX_ITER, diff);
        printf("Время выполнения: %f секунд\n", end_time - start_time);
    }


    for (i = 0; i < n_local + 2; i++) {
        free(f_old[i]);
        free(f_new[i]);
    }
    free(f_old);
    free(f_new);

    MPI_Finalize();
    return 0;
}

