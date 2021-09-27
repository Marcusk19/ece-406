#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
// #include <mpi.h>

#include "sort.h"
#include "graph.h"

struct Graph* countSortEdgesBySource (struct Graph* graph, int nbit){

    int i;
    int key;
    int pos;
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
    int radix = 8;

    // auxiliary arrays, allocated at the start up of the program
    int *vertex_count = (int*)malloc(graph->num_vertices*sizeof(int)); // needed for Counting Sort

    for(i = 0; i < graph->num_vertices; ++i) {
        vertex_count[i] = 0;
    }

    // count occurrence of key: id of a source vertex
    for(i = 0; i < graph->num_edges; ++i) {
        key = graph->sorted_edges_array[i].src;
        int temp = (key >> nbit*radix) & 0xffff;
        vertex_count[temp]++;
    }

    // transform to cumulative sum
    for(i = 1; i < graph->num_vertices; ++i) {
        vertex_count[i] += vertex_count[i - 1];
    }

    // fill-in the sorted array of edges
    for(i = graph->num_edges - 1; i >= 0; --i) {
        key = graph->sorted_edges_array[i].src;
        int temp = (key >> nbit*radix) & 0xffff;
        printf("temp:%d\n",temp);
        pos = vertex_count[temp] - 1;
        sorted_edges_array[pos] = graph->sorted_edges_array[i];
        vertex_count[temp]--;
    }



    free(vertex_count);
    free(graph->sorted_edges_array);

    graph->sorted_edges_array = sorted_edges_array;
    printEdgeArray(graph->sorted_edges_array, graph->num_edges);
    return graph;

}

// struct Graph* radixSortEdgesBySourceOpenMP (struct Graph* graph){

//         // printf("*** START Radix Sort Edges By Source OpenMP *** \n");

// }
// struct Graph* radixSortEdgesBySourceMPI (struct Graph* graph){

//         // printf("*** START Radix Sort Edges By Source MPI*** \n");

// }
// struct Graph* radixSortEdgesBySourceHybrid (struct Graph* graph){

//         // printf("*** START Radix Sort Edges By Source Hybrid*** \n");

// }
