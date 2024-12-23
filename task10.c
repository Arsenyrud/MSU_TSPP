#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define GRID_SIZE 128
#define NUM_ITERATIONS 100

#define INDEX(x, y, z) ((x) * (local_dim_y + 2) * (local_dim_x + 2) + (y) * (local_dim_x + 2) + (z))

int main(int argc, char **argv){
    int rank, size;
    int grid_dimensions[3], periods[3] = {0, 0, 0}, process_coords[3];
    MPI_Comm cart_comm;
    int local_dim_x, local_dim_y, local_dim_z;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    switch (size){
        case 1: grid_dimensions[0] = 1; grid_dimensions[1] = 1; grid_dimensions[2] = 1; break;
        case 2: grid_dimensions[0] = 1; grid_dimensions[1] = 1; grid_dimensions[2] = 2; break;
        case 4: grid_dimensions[0] = 1; grid_dimensions[1] = 2; grid_dimensions[2] = 2; break;
        case 6: grid_dimensions[0] = 1; grid_dimensions[1] = 2; grid_dimensions[2] = 3; break;
        case 8: grid_dimensions[0] = 2; grid_dimensions[1] = 2; grid_dimensions[2] = 2; break;
        case 12: grid_dimensions[0] = 2; grid_dimensions[1] = 3; grid_dimensions[2] = 2; break;
        case 16: grid_dimensions[0] = 2; grid_dimensions[1] = 4; grid_dimensions[2] = 2; break;
        default:
            if (rank == 0) {
                fprintf(stderr, "Unsupported number of processes.\n");
            }
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_Cart_create(MPI_COMM_WORLD, 3, grid_dimensions, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, rank, 3, process_coords);

    local_dim_x = GRID_SIZE / grid_dimensions[0];
    local_dim_y = GRID_SIZE / grid_dimensions[1];
    local_dim_z = GRID_SIZE / grid_dimensions[2];

    double *current_grid = (double*)malloc((local_dim_x + 2) * (local_dim_y + 2) * (local_dim_z + 2) * sizeof(double));
    double *updated_grid = (double*)malloc((local_dim_x + 2) * (local_dim_y + 2) * (local_dim_z + 2) * sizeof(double));
    if (!current_grid || !updated_grid) {
        fprintf(stderr, "Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    srand(time(NULL) + rank);
    for (int idx = 0; idx < (local_dim_x + 2) * (local_dim_y + 2) * (local_dim_z + 2); idx++){ 
        current_grid[idx] = (double)rand() / RAND_MAX;
        updated_grid[idx] = current_grid[idx];
    }

    MPI_Datatype yz_plane_type, xz_plane_type, xy_plane_type;
    MPI_Type_vector(local_dim_z, local_dim_y, (local_dim_y + 2)*(local_dim_x + 2), MPI_DOUBLE, &yz_plane_type);
    MPI_Type_create_resized(yz_plane_type, 0, local_dim_y * sizeof(double), &yz_plane_type);
    MPI_Type_commit(&yz_plane_type);

    MPI_Type_vector(local_dim_z, local_dim_x, local_dim_x + 2, MPI_DOUBLE, &xz_plane_type);
    MPI_Type_create_resized(xz_plane_type, 0, local_dim_x * sizeof(double), &xz_plane_type);
    MPI_Type_commit(&xz_plane_type);

    MPI_Type_contiguous(local_dim_x * local_dim_y, MPI_DOUBLE, &xy_plane_type);
    MPI_Type_create_resized(xy_plane_type, 0, local_dim_x * local_dim_y * sizeof(double), &xy_plane_type);
    MPI_Type_commit(&xy_plane_type);

    int north_neighbor, south_neighbor, east_neighbor, west_neighbor, up_neighbor, down_neighbor;
    MPI_Cart_shift(cart_comm, 0, 1, &west_neighbor, &east_neighbor);
    MPI_Cart_shift(cart_comm, 1, 1, &south_neighbor, &north_neighbor);
    MPI_Cart_shift(cart_comm, 2, 1, &down_neighbor, &up_neighbor);

    double start_time = MPI_Wtime();

    for (int iter = 0; iter < NUM_ITERATIONS; iter++){
        MPI_Request requests[12];
        int req_count = 0;

        if (west_neighbor != MPI_PROC_NULL){
            MPI_Isend(&current_grid[INDEX(1,1,1)], 1, yz_plane_type, west_neighbor, 0, cart_comm, &requests[req_count++]);
            MPI_Irecv(&current_grid[INDEX(1,1,0)], 1, yz_plane_type, west_neighbor, 1, cart_comm, &requests[req_count++]);
        }
        if (east_neighbor != MPI_PROC_NULL){
            MPI_Isend(&current_grid[INDEX(1,1,local_dim_x)], 1, yz_plane_type, east_neighbor, 1, cart_comm, &requests[req_count++]);
            MPI_Irecv(&current_grid[INDEX(1,1,local_dim_x+1)], 1, yz_plane_type, east_neighbor, 0, cart_comm, &requests[req_count++]);
        }

        if (south_neighbor != MPI_PROC_NULL){
            MPI_Isend(&current_grid[INDEX(1,1,1)], 1, xz_plane_type, south_neighbor, 2, cart_comm, &requests[req_count++]);
            MPI_Irecv(&current_grid[INDEX(1,0,1)], 1, xz_plane_type, south_neighbor, 3, cart_comm, &requests[req_count++]);
        }
        if (north_neighbor != MPI_PROC_NULL){
            MPI_Isend(&current_grid[INDEX(1,local_dim_y,1)], 1, xz_plane_type, north_neighbor, 3, cart_comm, &requests[req_count++]);
            MPI_Irecv(&current_grid[INDEX(1,local_dim_y+1,1)], 1, xz_plane_type, north_neighbor, 2, cart_comm, &requests[req_count++]);
        }

        if (down_neighbor != MPI_PROC_NULL){
            MPI_Isend(&current_grid[INDEX(1,1,1)], 1, xy_plane_type, down_neighbor, 4, cart_comm, &requests[req_count++]);
            MPI_Irecv(&current_grid[INDEX(0,1,1)], 1, xy_plane_type, down_neighbor, 5, cart_comm, &requests[req_count++]);
        }
        if (up_neighbor != MPI_PROC_NULL){
            MPI_Isend(&current_grid[INDEX(local_dim_z,1,1)], 1, xy_plane_type, up_neighbor, 5, cart_comm, &requests[req_count++]);
            MPI_Irecv(&current_grid[INDEX(local_dim_z+1,1,1)], 1, xy_plane_type, up_neighbor, 4, cart_comm, &requests[req_count++]);
        }

        MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);

        for (int z = 1; z <= local_dim_z; z++){
            for (int y = 1; y <= local_dim_y; y++){
                for (int x = 1; x <= local_dim_x; x++){
                    updated_grid[INDEX(z,y,x)] = ( 
                        current_grid[INDEX(z-1,y,x)] + current_grid[INDEX(z+1,y,x)] +
                        current_grid[INDEX(z,y-1,x)] + current_grid[INDEX(z,y+1,x)] +
                        current_grid[INDEX(z,y,x-1)] + current_grid[INDEX(z,y,x+1)]
                    ) / 6.0;
                }
            }
        }

        if (iter != NUM_ITERATIONS -1){
            for (int z = 1; z <= local_dim_z; z++){
                for (int y = 1; y <= local_dim_y; y++){
                    for (int x = 1; x <= local_dim_x; x++){
                        current_grid[INDEX(z,y,x)] = updated_grid[INDEX(z,y,x)];
                    }
                }
            }
        }
    }

    double end_time = MPI_Wtime();

    double local_difference = 0.0;
    for (int z = 1; z <= local_dim_z; z++){
        for (int y = 1; y <= local_dim_y; y++){
            for (int x = 1; x <= local_dim_x; x++){
                double diff = updated_grid[INDEX(z,y,x)] - current_grid[INDEX(z,y,x)];
                local_difference += diff * diff;
            }
        }
    }

    double global_difference;
    MPI_Reduce(&local_difference, &global_difference, 1, MPI_DOUBLE, MPI_SUM, 0, cart_comm);

    if (rank == 0){
        printf("Норма: %lf\n", sqrt(global_difference / (GRID_SIZE * GRID_SIZE * GRID_SIZE)));
        printf("Время выполнения: %lf секунд\n", end_time - start_time);
    }

    MPI_Type_free(&yz_plane_type);
    MPI_Type_free(&xz_plane_type);
    MPI_Type_free(&xy_plane_type);

    free(current_grid);
    free(updated_grid);

    MPI_Finalize();
    return 0;
}

