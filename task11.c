#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N 10000
#define PRINT_COUNT 5

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int dims[2];
    switch(size) {
        case 1:  dims[0] = 1; dims[1] = 1; break;
        case 2:  dims[0] = 2; dims[1] = 1; break;
        case 4:  dims[0] = 2; dims[1] = 2; break;
        case 6:  dims[0] = 3; dims[1] = 2; break;
        case 8:  dims[0] = 4; dims[1] = 2; break;
        case 12: dims[0] = 4; dims[1] = 3; break;
        case 16: dims[0] = 4; dims[1] = 4; break;
        default: dims[0] = size; dims[1] = 1;
    }

    int row_coord = rank / dims[1];
    int col_coord = rank % dims[1];
    int Ny = N / dims[0];
    int Nx = N / dims[1];
    double *A = (double *)malloc(Ny * Nx * sizeof(double));
    double *c = (double *)calloc(Ny, sizeof(double));
    srand(time(NULL) + rank);
    for (int i = 0; i < Ny * Nx; i++) {
        A[i] = (double)rand() / RAND_MAX;
    }

    printf("Process %d/%d: row_coord=%d col_coord=%d, Ny=%d Nx=%d\n",
           rank, size, row_coord, col_coord, Ny, Nx);
    printf("Process %d: A[0..%d] = ", rank, PRINT_COUNT - 1);
    for (int i = 0; i < Ny * Nx && i < PRINT_COUNT; i++) {
        printf("%.3f ", A[i]);
    }
    printf("\n");
    fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);

    double *b = NULL;
    MPI_Win win_b;
    MPI_Win_allocate(N * sizeof(double), sizeof(double),
                     MPI_INFO_NULL, MPI_COMM_WORLD, &b, &win_b);

    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            b[i] = (double)rand() / RAND_MAX;
        }
    }

    MPI_Win_fence(0, win_b);

    if (rank == 0) {
        printf("Global b (first %d elements): ", PRINT_COUNT);
        for (int i = 0; i < N && i < PRINT_COUNT; i++) {
            printf("%.3f ", b[i]);
        }
        printf("\n");
        fflush(stdout);
    }

    double *b_local = (double *)malloc(Nx * sizeof(double));
    MPI_Get(b_local, Nx, MPI_DOUBLE, 0, col_coord * Nx, Nx, MPI_DOUBLE, win_b);
    MPI_Win_fence(0, win_b);

    printf("Process %d: b_local[0..%d] = ", rank, PRINT_COUNT - 1);
    for (int i = 0; i < Nx && i < PRINT_COUNT; i++) {
        printf("%.3f ", b_local[i]);
    }
    printf("\n");
    fflush(stdout);

    MPI_Barrier(MPI_COMM_WORLD);

    double t_start = MPI_Wtime();
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            c[i] += A[i * Nx + j] * b_local[j];
        }
    }
    double t_calc = MPI_Wtime() - t_start;

    printf("Process %d: c[0..%d] = ", rank, PRINT_COUNT - 1);
    for (int i = 0; i < Ny && i < PRINT_COUNT; i++) {
        printf("%.3f ", c[i]);
    }
    printf("\n");
    fflush(stdout);

    double *c_global = NULL;
    MPI_Win win_c;
    MPI_Win_allocate(N * sizeof(double), sizeof(double),
                     MPI_INFO_NULL, MPI_COMM_WORLD, &c_global, &win_c);

    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            c_global[i] = 0.0;
        }
    }

    MPI_Win_fence(0, win_c);
    MPI_Win_lock_all(0, win_c);
    MPI_Accumulate(c, Ny, MPI_DOUBLE, 0, row_coord * Ny, Ny, MPI_DOUBLE, MPI_SUM, win_c);
    MPI_Win_unlock_all(win_c);
    MPI_Win_fence(0, win_c);

    double t_elapsed = MPI_Wtime() - t_start;

    if (rank == 0) {
        printf("==== Final results on rank 0 ====\n");
        printf("c_global[0..%d] = ", PRINT_COUNT - 1);
        for (int i = 0; i < N && i < PRINT_COUNT; i++) {
            printf("%.3f ", c_global[i]);
        }
        printf("\n");
        printf("Time (total) = %lf\n", t_elapsed);
        printf("Time (calc only) = %lf\n", t_calc);
    }

    MPI_Win_free(&win_c);
    MPI_Win_free(&win_b);
    free(b_local);
    free(A);
    free(c);

    MPI_Finalize();
    return 0;
}

