#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

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

struct Graph *countSortEdgesBySource (struct Graph *graph, int nbit)
{   
    printf("count sort for nbit %d\n", nbit);
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
    // auxillary arrays, allocated at the start up of the program
    u_int32_t *vertex_count = (u_int32_t *)malloc(graph->num_vertices * sizeof(u_int32_t)); //neede for coutning sort
    #pragma omp parallel
    {
        int i;
        int key;
        int pos;
        int P = omp_get_num_threads();
        int tid = omp_get_thread_num(); // used for parallization
        int offset_start, offset_end; // used for parallization
        int vert_start, vert_end;
        int radix = 2; // used for radix sort
        if(tid == 0){
           vertex_count = (u_int32_t*)malloc(P*graph->num_vertices*sizeof(u_int32_t));
        }
        //printf("malloc:%d\n", (int)P*graph->num_edges*sizeof(u_int32_t));
        
        #pragma omp barrier

        offset_start = tid * (graph->num_edges/P);
        vert_start = tid * (graph->num_vertices/P);
        if(tid == (P-1)){
            offset_end = offset_start + (graph->num_edges/P) + (graph->num_edges%P); // edge case for last group of indices
            vert_end = vert_start + (graph->num_vertices/P) + (graph->num_vertices%P);
        }else{
            offset_end = offset_start + (graph->num_edges/P);
            vert_end = vert_start + (graph->num_vertices/P);
        }
        printf("num_vertices:%d\n", graph->num_vertices);
        printf("num_edges:%d\n", graph->num_edges);
        printf("%d to %d\n", offset_start, offset_end);
        #pragma omp barrier

        // auxiliary arrays, allocated at the start up of the program
        //int *vertex_count = (int *)malloc(graph->num_vertices * sizeof(int)); // needed for Counting Sort

        printf("Started array initialization\n");
        for(i = vert_start; i < vert_end; ++i)
        {
            //printf("%d\n", i);
            vertex_count[i] = 0; 
            //printf("vertex_count[%d]=%d\n",i,vertex_count[i]);       
        }
        //printf("Finished array initialization\n");
        #pragma omp barrier

        // count occurrence of key: id of a source vertex
        printf("Started count occurence\n");
        for(i = offset_start; i < offset_end; ++i)
        {
            key = graph->sorted_edges_array[i].src;
            key = (key >> (radix * nbit)) & 0xff;
            //printf("nbit = %d key = %d\n", nbit, key);
            vertex_count[key]++;
            // printf("vertex_count[%d]=%d\n", key,vertex_count[key]);
        }
        //printf("Finished count occurence\n");
        #pragma omp barrier

        // transform to cumulative sum
        printf("Started cumulative sum\n");
        for(i = vert_start+1; i < vert_end; ++i)  
        {
            vertex_count[i] += vertex_count[i - 1];
            //printf("vertex_count[%d]=%d\n",i,vertex_count[i]);
        }
        //printf("Finished cumulative sum\n");
        #pragma omp barrier

        // fill-in the sorted array of edges
        printf("Started fill-in\n");
        for(i = offset_end - 1; i >= offset_start; --i)
        {
            key = graph->sorted_edges_array[i].src;
            //printf("1\n");
            key = (key >> (radix * nbit)) & 0xff; // find 8-bit radix
            //printf("key:%d\n", key);
            pos = vertex_count[key] - 1;
            //printf("pos:%d\n", pos);
            sorted_edges_array[pos] = graph->sorted_edges_array[i];
            //printf("4\n");
            vertex_count[key]--;
        }  
        printf("Finished fill-in\n");
        #pragma omp barrier
    } //end parallel region

    free(vertex_count);
    free(graph->sorted_edges_array);
    graph->sorted_edges_array = sorted_edges_array;
    return graph; 
}

#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{
    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
    int maxSrc = graph->sorted_edges_array[0].src;
    for(int i = 1; i < graph->num_edges; i++){
        if(graph->sorted_edges_array[i].src > maxSrc) maxSrc = graph->sorted_edges_array[i].src;
    }
    printf("maxSrc = %d\n", maxSrc);
    int n = sizeof(maxSrc);
    for(int i = 1; i <= n; i++){
        graph = countSortEdgesBySource(graph, i);
        //printEdgeArray(graph->sorted_edges_array, graph->num_edges);
    }
    //printEdgeArray(graph->sorted_edges_array, graph->num_edges);
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
