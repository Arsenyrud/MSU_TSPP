#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <time.h>
#include <math.h> 


void matmul_seq(double *A, double *B, double *C, int N);
void matmul_avx(double *A, double *B, double *C, int N);

int main() {
    int N = 2048;
    double *A, *B, *C_seq, *C_avx;
    posix_memalign((void**)&A, 32, N * N * sizeof(double));
    posix_memalign((void**)&B, 32, N * N * sizeof(double));
    posix_memalign((void**)&C_seq, 32, N * N * sizeof(double));
    posix_memalign((void**)&C_avx, 32, N * N * sizeof(double));

    for (int i = 0; i < N * N; i++) {
        A[i] = (double)rand() / RAND_MAX;
        B[i] = (double)rand() / RAND_MAX;
        C_seq[i] = 0.0;
        C_avx[i] = 0.0;
    }

    struct timespec start, end;

    // Последовательная версия
    clock_gettime(CLOCK_MONOTONIC, &start);
    matmul_seq(A, B, C_seq, N);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time_seq = (end.tv_sec - start.tv_sec) + 
                      (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Sequential time: %f seconds\n", time_seq);

    // Векторизованная версия
    clock_gettime(CLOCK_MONOTONIC, &start);
    matmul_avx(A, B, C_avx, N);
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time_avx = (end.tv_sec - start.tv_sec) + 
                      (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("AVX time: %f seconds\n", time_avx);

    double max_diff = 0.0;
    for (int i = 0; i < N * N; i++) {
        double diff = fabs(C_seq[i] - C_avx[i]);
        if (diff > max_diff) {
            max_diff = diff;
        }
    }
    printf("Max difference: %e\n", max_diff);

    free(A);
    free(B);
    free(C_seq);
    free(C_avx);

    return 0;
}

void matmul_seq(double *A, double *B, double *C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k < N; k++) {
                sum += A[i*N + k] * B[k*N + j];
            }
            C[i*N + j] = sum;
        }
    }
}

void matmul_avx(double *A, double *B, double *C, int N) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j += 4) {
            __m256d c_vec = _mm256_setzero_pd();
            for (int k = 0; k < N; k++) {
                __m256d a_vec = _mm256_broadcast_sd(&A[i*N + k]);
                __m256d b_vec = _mm256_load_pd(&B[k*N + j]);
                c_vec = _mm256_fmadd_pd(a_vec, b_vec, c_vec);
            }
            _mm256_store_pd(&C[i*N + j], c_vec);
        }
    }
}