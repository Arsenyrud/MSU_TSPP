#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define WIDTH  100
#define HEIGHT 100
#define K 10

void initialize_grid(int *grid, int width, int height) {
    for (int i = 0; i < width * height; i++) {
        grid[i] = rand() % 2;
    }
}

int count_neighbors(int *grid, int x, int y, int width, int height) {
    int count = 0;
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            if (dx != 0 || dy != 0) {
                int nx = (x + dx + width) % width;
                int ny = (y + dy + height) % height;
                count += grid[ny * width + nx];
            }
        }
    }
    return count;
}

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);                
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
    MPI_Comm_size(MPI_COMM_WORLD, &size);  

    int local_height = HEIGHT / size;
    int remainder = HEIGHT % size;
    if (rank < remainder) {
        local_height += 1;
    }

    int global_y_start;
    if (rank < remainder) {
        global_y_start = rank * local_height;
    } else {
        global_y_start = rank * local_height + remainder;
    }

    int *local_grid = (int *)malloc(WIDTH * (local_height + 2) * sizeof(int));
    int *new_local_grid = (int *)malloc(WIDTH * (local_height + 2) * sizeof(int));

    int *full_grid = NULL;
    int *sendcounts = NULL;
    int *displs = NULL;
    if (rank == 0) {
        full_grid = (int *)malloc(WIDTH * HEIGHT * sizeof(int));
        initialize_grid(full_grid, WIDTH, HEIGHT);

        sendcounts = (int *)malloc(size * sizeof(int));
        displs = (int *)malloc(size * sizeof(int));

        int offset = 0;
        for (int i = 0; i < size; i++) {
            int rows = HEIGHT / size;
            if (i < remainder) {
                rows += 1;
            }
            sendcounts[i] = rows * WIDTH;
            displs[i] = offset;
            offset += rows * WIDTH;
        }
    }

    int recvcount = local_height * WIDTH;
    MPI_Scatterv(full_grid, sendcounts, displs, MPI_INT,
                 &local_grid[WIDTH], recvcount, MPI_INT,
                 0, MPI_COMM_WORLD);

    if (rank == 0) {
        free(full_grid);
        free(sendcounts);
        free(displs);
    }

    int iteration = 0;
    int stable = 0;
    int prev_total_live_cells = -1;

    double start_time = MPI_Wtime();

    while (!stable) {
        iteration++;

        MPI_Request requests[4];
        int up_rank = (rank - 1 + size) % size;
        int down_rank = (rank + 1) % size;

        MPI_Isend(&local_grid[WIDTH], WIDTH, MPI_INT, up_rank, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Irecv(&local_grid[0], WIDTH, MPI_INT, up_rank, 1, MPI_COMM_WORLD, &requests[1]);

        MPI_Isend(&local_grid[local_height * WIDTH], WIDTH, MPI_INT, down_rank, 1, MPI_COMM_WORLD, &requests[2]);
        MPI_Irecv(&local_grid[(local_height + 1) * WIDTH], WIDTH, MPI_INT, down_rank, 0, MPI_COMM_WORLD, &requests[3]);

        MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

        int local_live_cells = 0;
        for (int y = 1; y <= local_height; y++) {
            for (int x = 0; x < WIDTH; x++) {
                int neighbors = count_neighbors(local_grid, x, y, WIDTH, local_height + 2);
                int idx = y * WIDTH + x;
                if (local_grid[idx] == 1) {
                    if (neighbors < 2 || neighbors > 3) {
                        new_local_grid[idx] = 0;
                    } else {
                        new_local_grid[idx] = 1;
                        local_live_cells++;
                    }
                } else {
                    if (neighbors == 3) {
                        new_local_grid[idx] = 1;
                        local_live_cells++;
                    } else {
                        new_local_grid[idx] = 0;
                    }
                }
            }
        }

        memcpy(&local_grid[WIDTH], &new_local_grid[WIDTH], WIDTH * local_height * sizeof(int));

        int total_live_cells;
        MPI_Allreduce(&local_live_cells, &total_live_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if (iteration >= K) {
            stable = (prev_total_live_cells == total_live_cells);
        }

        int global_stable;
        MPI_Allreduce(&stable, &global_stable, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        stable = global_stable;

        prev_total_live_cells = total_live_cells;

        if (stable) {
            break;
        }
    }

    double end_time = MPI_Wtime();
    double total_time = end_time - start_time;

    if (rank == 0) {
        printf("Игра остановлена на итерации %d. Общее число живых клеток: %d\n", iteration, prev_total_live_cells);
        printf("Время выполнения: %f секунд\n", total_time);
    }

    free(local_grid);
    free(new_local_grid);
    MPI_Finalize();
    return 0;
}
