#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "papi.h"

#define N_TESTS 5

typedef struct
{
    int row_count;
    unsigned int col_count;
    unsigned int *row_ptr;
    int *col_ids;
    double *vals;
} CSR_graph;

int read_csr_graph(const char *filename, CSR_graph *graph)
{
    FILE *graph_file = fopen(filename, "rb");

    int fd = fread(&graph->row_count, sizeof(int), 1, graph_file);

    fd = fread(&graph->col_count, sizeof(unsigned int), 1, graph_file);

    printf("Чтение графа из файла: %s\n", filename);
    printf("Количество вершин: %d, количество рёбер: %u\n", graph->row_count, graph->col_count);

    graph->row_ptr = (unsigned int *)malloc((graph->row_count + 1) * sizeof(unsigned int));
    graph->col_ids = (int *)malloc(graph->col_count * sizeof(int));
    graph->vals = (double *)malloc(graph->col_count * sizeof(double));

    if (fread(graph->row_ptr, sizeof(unsigned int), graph->row_count + 1, graph_file) != (size_t)(graph->row_count + 1))
    {
        fclose(graph_file);
        return -1;
    }

    if (fread(graph->col_ids, sizeof(int), graph->col_count, graph_file) != graph->col_count)
    {
        fclose(graph_file);
        return -1;
    }

    if (fread(graph->vals, sizeof(double), graph->col_count, graph_file) != graph->col_count)
    {
        fclose(graph_file);
        return -1;
    }

    fclose(graph_file);
    return 0;
}

void print_vertex(const CSR_graph *graph, int idx)
{
    printf("Вершина %d инцидентные рёбра:\n", idx);
    for (unsigned int col = graph->row_ptr[idx]; col < graph->row_ptr[idx + 1]; col++)
    {
        printf("  -> %d с весом %.2lf\n", graph->col_ids[col], graph->vals[col]);
    }
    printf("\n");
}

void reset_graph(CSR_graph *graph)
{
    graph->row_count = 0;
    graph->col_count = 0;
    free(graph->row_ptr);
    free(graph->col_ids);
    free(graph->vals);
    graph->row_ptr = NULL;
    graph->col_ids = NULL;
    graph->vals = NULL;
}

int algorithm1(const CSR_graph *graph)
{
    int max_vertex = -1;
    double max_total_weight = -1.0;

    for (int v = 0; v < graph->row_count; v++)
    {
        double total_weight = 0.0;
        unsigned int start_edge = graph->row_ptr[v];
        unsigned int end_edge = graph->row_ptr[v + 1];

        for (unsigned int e = start_edge; e < end_edge; e++)
        {
            int dest_v = graph->col_ids[e];
            if (dest_v % 2 == 0)
            {
                total_weight += graph->vals[e];
            }
        }

        if (total_weight > max_total_weight)
        {
            max_total_weight = total_weight;
            max_vertex = v;
        }
    }

    printf("Алгоритм 1: Вершина %d имеет максимальный суммарный вес %.6lf\n", max_vertex, max_total_weight);
    return max_vertex;
}

int algorithm2(const CSR_graph *graph)
{
    int max_vertex = -1;
    double max_rank = -1.0;

    for (int v = 0; v < graph->row_count; v++)
    {
        double rank = 0.0;
        unsigned int start_edge_v = graph->row_ptr[v];
        unsigned int end_edge_v = graph->row_ptr[v + 1];

        for (unsigned int e_v = start_edge_v; e_v < end_edge_v; e_v++)
        {
            double w_edge_i = graph->vals[e_v];
            int dest_v_i = graph->col_ids[e_v];

            double W_dest_v_i = 0.0;
            unsigned int start_edge_dest = graph->row_ptr[dest_v_i];
            unsigned int end_edge_dest = graph->row_ptr[dest_v_i + 1];

            for (unsigned int e_dest = start_edge_dest; e_dest < end_edge_dest; e_dest++)
            {
                double w_edge_j = graph->vals[e_dest];
                int dest_v_j = graph->col_ids[e_dest];

                int N_inc_edges_dest_vj = graph->row_ptr[dest_v_j + 1] - graph->row_ptr[dest_v_j];
                W_dest_v_i += w_edge_j * N_inc_edges_dest_vj;
            }

            rank += w_edge_i * W_dest_v_i;
        }

        if (rank > max_rank)
        {
            max_rank = rank;
            max_vertex = v;
        }
    }

    printf("Алгоритм 2: Вершина %d имеет наивысший ранг %.6lf\n", max_vertex, max_rank);
    return max_vertex;
}

void init_papi(int *EventSet)
{
    int retval;
    *EventSet = PAPI_NULL;

    retval = PAPI_library_init(PAPI_VER_CURRENT);
    retval = PAPI_create_eventset(EventSet);
    retval = PAPI_add_event(*EventSet, PAPI_L1_TCM);
    retval = PAPI_add_event(*EventSet, PAPI_L2_TCM);
    retval = PAPI_add_event(*EventSet, PAPI_TOT_INS);
}

void measure_algorithm(int (*algorithm)(const CSR_graph *), const CSR_graph *graph, const char *algo_name)
{
    int EventSet = PAPI_NULL;
    long long values[3];
    int retval;

    init_papi(&EventSet);
    retval = PAPI_start(EventSet);

    int result = algorithm(graph);

    retval = PAPI_stop(EventSet, values);

    printf("\nПоказатели PAPI %s:\n", algo_name);
    printf("PAPI_L1_TCM: %lld\n", values[0]);
    printf("PAPI_L2_TCM: %lld\n", values[1]);
    printf("PAPI_TOT_INS: %lld\n", values[2]);

    PAPI_cleanup_eventset(EventSet);
    PAPI_destroy_eventset(&EventSet);
    PAPI_shutdown();
}

int main()
{
    const char *filenames[N_TESTS] = {"synt", "road_graph", "stanford", "youtube", "syn_rmat"};

    for (int n_test = 0; n_test < N_TESTS; n_test++)
    {
        printf("\n===============================\n");
        printf("Тест #%d: %s\n", n_test + 1, filenames[n_test]);
        printf("===============================\n");

        CSR_graph graph;
        read_csr_graph(filenames[n_test], &graph);

        measure_algorithm(algorithm1, &graph, "Алгоритм 1");
        measure_algorithm(algorithm2, &graph, "Алгоритм 2");

        reset_graph(&graph);
    }

    return 0;
}