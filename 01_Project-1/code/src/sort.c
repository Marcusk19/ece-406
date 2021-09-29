#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>

#ifdef OPENMP_HARNESS
#include <omp.h>
#endif

#ifdef MPI_HARNESS
#include <mpi.h>
#endif

#ifdef HYBRID_HARNESS
#include <omp.h>
#include <mpi.h>
#endif

#include "sort.h"
#include "graph.h"
#include "edgelist.h"

struct Graph *countSortEdgesBySource (struct Graph *graph, int radix)
{   
    printf("count sort for radix %d\n", radix);
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
    int i;
    int key;
    int pos;
    int P;
    int nbit = 8; // used for radix sort
    int base = 0;

    // auxillary arrays, allocated at the start up of the program
    int *vertex_count; // needed for coutning sort
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int offset_start, offset_end;
        printf("tid: %d\n", tid);
        if(tid == 0){
            P = omp_get_num_threads();
            vertex_count = (int*)malloc(P*graph->num_vertices*sizeof(int));
        }
        #pragma omp barrier
        offset_start = tid * (graph->num_edges/P);
        if(tid == (P-1)){
            offset_end = offset_start + (graph->num_edges/P) + (graph->num_edges%P); // edge case for last group of indices
        }else{
            offset_end = offset_start + (graph->num_edges/P);
        }
        int tid_start = tid * (graph->num_vertices);
        int tid_end = tid * (graph->num_vertices) + graph->num_vertices;
        if(tid == 0){
            printf("num edges: %d\n", graph->num_edges);
            printf("num verts: %d\n", graph->num_vertices);
        }
        printf("edges %d to %d\n", offset_start, offset_end);
        printf("vertices %d to %d\n", tid_start, tid_end);

        for(i = tid_start; i < tid_end; i++){
            vertex_count[i] = 0;
        }
 
        printf("tid: %d done\n", tid);
        // count occurrence of key: id of a source vertex
        printf("Started count occurence\n");
        for(int i = offset_start; i < offset_end; ++i)
        {
            key = graph->sorted_edges_array[i].src;
            int temp = (key>>(radix * nbit)) & 0xff;
            vertex_count[tid * graph->num_vertices + temp]++;
        }
        printf("vertex_count[%d]: %d\n", tid_start, vertex_count[tid_start]);
        #pragma omp barrier
        
        if(tid == 0){
            for(int i = 0; i < graph->num_vertices; i++){
                    for(int t = 0; t < P; ++t){
                            base += vertex_count[t*graph->num_vertices + i];
                            vertex_count[t*graph->num_vertices + i] = base;
                    } 
            }
        }
                
        #pragma omp barrier
        printf("vertex_count[%d]: %d\n", tid_end-1, vertex_count[tid_end-1]);
                                                                            
        
        for(int i = offset_end - 1; i >= offset_start; i--){
            key = graph->sorted_edges_array[i].src;
            int temp = (key>>(radix * nbit)) & 0xff;
            pos = vertex_count[tid*graph->num_vertices+temp]-1;
            sorted_edges_array[pos] = graph->sorted_edges_array[i];
            vertex_count[tid*graph->num_vertices+temp]--;
        }
    } //end parallel region
    printf("end parallel region\n");
    free(vertex_count);
    free(graph->sorted_edges_array);
    graph->sorted_edges_array = sorted_edges_array;
    //printEdgeArray(graph->sorted_edges_array, graph->num_edges);
    return graph; 
}

#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{
    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
    int n = sizeof(graph->num_vertices);
    int P;
    printf("n:%d \n", n);
    #pragma omp parallel
    { // start parallel region
        int offset_start;
        int offset_end;
        int tid = omp_get_thread_num();
        if(tid == 0){
            P = omp_get_num_threads();
        }

        offset_start = tid * (n/P);
    } // end parallel region
    for(int i = 1; i <= 4; i++){
        graph = countSortEdgesBySource(graph, i);
    }
    printEdgeArray(graph->sorted_edges_array, 100);
    return graph;
}
#endif

#ifdef MPI_HARNESS
struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source MPI*** \n");
    return graph;
}
#endif

#ifdef HYBRID_HARNESS
struct Graph *radixSortEdgesBySourceHybrid (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source Hybrid*** \n");
    return graph;
}
#endif
