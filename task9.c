#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

void initialize_matrix(int *mat, int N, int value) {
    for (int i = 0; i < N * N; i++) {
        mat[i] = value;
    }
}

void multiply_submatrices(int *A, int *B, int *C, int block_size) {
    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < block_size; j++) {
            for (int k = 0; k < block_size; k++) {
                C[i * block_size + j] += A[i * block_size + k] * B[k * block_size + j];
            }
        }
    }
}

int verify_result(int *C, int N, int expected) {
    for (int i = 0; i < N * N; i++) {
        if (C[i] != expected) {
            return 0;
        }
    }
    return 1;
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Comm cart_comm, row_comm, col_comm;
    int dims[2], periods[2] = {0, 0}, coords[2];
    int reorder = 1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 1000; // Размер матрицы
    int b = 100;  // Размер блока

    int q = (int)sqrt((double)size);
    if (q * q != size) {
        if (rank == 0) {
            printf("Number of processes must be a perfect square.\n");
        }
        MPI_Finalize();
        exit(0);
    }

    dims[0] = dims[1] = q; 
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    MPI_Comm_split(cart_comm, coords[0], coords[1], &row_comm);
    MPI_Comm_split(cart_comm, coords[1], coords[0], &col_comm);

    if (N % (q * b) != 0) {
        if (rank == 0) {
            printf("Matrix size N must be divisible by sqrt(P) * b.\n");
        }
        MPI_Finalize();
        exit(0);
    }

    int block_size = N / q;

    int *A_block = (int *)malloc(block_size * block_size * sizeof(int));
    int *B_block = (int *)malloc(block_size * block_size * sizeof(int));
    int *C_block = (int *)malloc(block_size * block_size * sizeof(int));
    memset(C_block, 0, block_size * block_size * sizeof(int));

    int *A = NULL;
    int *B = NULL;
    int *C = NULL;

    if (rank == 0) {
        A = (int *)malloc(N * N * sizeof(int));
        B = (int *)malloc(N * N * sizeof(int));
        C = (int *)malloc(N * N * sizeof(int));
        initialize_matrix(A, N, 1);  
        initialize_matrix(B, N, 1);  
    }


    MPI_Datatype block_type, resized_block_type;
    MPI_Type_create_subarray(2, (int[]){N, N}, (int[]){block_size, block_size},
                             (int[]){coords[0] * block_size, coords[1] * block_size}, MPI_ORDER_C, MPI_INT, &block_type);
    MPI_Type_create_resized(block_type, 0, block_size * sizeof(int), &resized_block_type);
    MPI_Type_commit(&resized_block_type);

    int *recvcounts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        recvcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));
        for (int i = 0; i < size; i++) {
            recvcounts[i] = 1;
        }
        for (int i = 0; i < size; i++) {
            displs[i] = i;
        }
    }


    MPI_Scatter(A, 1, resized_block_type, A_block, block_size * block_size, MPI_INT, 0, cart_comm);
    MPI_Scatter(B, 1, resized_block_type, B_block, block_size * block_size, MPI_INT, 0, cart_comm);

    MPI_Type_free(&block_type);
    MPI_Type_free(&resized_block_type);
    if (rank == 0) {
        free(recvcounts);
        free(displs);
    }


    int *A_temp = (int *)malloc(block_size * block_size * sizeof(int));
    int *B_temp = (int *)malloc(block_size * block_size * sizeof(int));

    double start_time = MPI_Wtime();


    for (int k = 0; k < q; k++) {
        if (coords[1] == k) {
            memcpy(A_temp, A_block, block_size * block_size * sizeof(int));
        }
        MPI_Bcast(A_temp, block_size * block_size, MPI_INT, k, row_comm);

        if (coords[0] == k) {
            memcpy(B_temp, B_block, block_size * block_size * sizeof(int));
        }
        MPI_Bcast(B_temp, block_size * block_size, MPI_INT, k, col_comm);

        multiply_submatrices(A_temp, B_temp, C_block, block_size);
    }

    double end_time = MPI_Wtime();

    free(A_temp);
    free(B_temp);

    MPI_Gather(C_block, block_size * block_size, MPI_INT, C, block_size * block_size, MPI_INT, 0, cart_comm);

    if (rank == 0) {
        printf("\nMatrix C (Result):\n");
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                printf("%d ", C[i * N + j]);
            }
            printf("\n");
        }

        int correct = verify_result(C, N, N);
        if (correct) {
            printf("\nResult is correct.\n");
        } else {
            printf("\nResult is incorrect.\n");
        }

        printf("\nExecution time: %f seconds.\n", end_time - start_time);

        free(A);
        free(B);
        free(C);
    }

    free(A_block);
    free(B_block);
    free(C_block);

    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Comm_free(&cart_comm);

    MPI_Finalize();
    return 0;
}

