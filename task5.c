#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>


void print_usage(char *prog_name) {
    printf("Usage: %s a b x p N P\n", prog_name);
    printf("Where:\n");
    printf("  a - left boundary (integer)\n");
    printf("  b - right boundary (integer), must satisfy a < x < b\n");
    printf("  x - starting point (integer)\n");
    printf("  p - probability to move right (float between 0 and 1)\n");
    printf("  N - number of particles (positive integer)\n");
    printf("  P - number of threads (positive integer)\n");
}

int main(int argc, char *argv[]) {
    if (argc != 7) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    int a = atoi(argv[1]);
    int b = atoi(argv[2]);
    int x = atoi(argv[3]);
    double p = atof(argv[4]);
    long long N = atoll(argv[5]);
    int P = atoi(argv[6]);

    if (!(a < x && x < b)) {
        return EXIT_FAILURE;
    }
    if (!(p >= 0.0 && p <= 1.0)) {
        return EXIT_FAILURE;
    }
    if (N <= 0) {
        return EXIT_FAILURE;
    }
    if (P <= 0) {
        return EXIT_FAILURE;
    }

    omp_set_num_threads(P);

    long long count_b = 0;
    long long total_steps = 0;

    double loop_start = omp_get_wtime();

    #pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        unsigned int seed = (unsigned int)(time(NULL)) + thread_num;

        #pragma omp for reduction(+: count_b, total_steps)
        for (long long i = 0; i < N; i++) {
            int pos = x;
            long long steps = 0;
            while (pos != a && pos != b) {
                double rand_val = (double)rand_r(&seed) / RAND_MAX;
                if (rand_val < p) {
                    pos += 1;
                } else {
                    pos -= 1;
                }
                steps++;
            }

            if (pos == b) {
                count_b += 1;
            }
            total_steps += steps;
        }
    }

    double loop_end = omp_get_wtime();
    double loop_time = loop_end - loop_start;

    double probability_b = ((double)count_b) / N;

    double average_steps = ((double)total_steps) / N;

    printf("Probability of reaching %d: %.6f\n", b, probability_b);
    printf("Average number of steps: %.2f\n", average_steps);
    printf("Time of the main loop: %.6f seconds\n", loop_time);

    return EXIT_SUCCESS;
}
