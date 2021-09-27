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

struct Graph *countSortEdgesBySource (struct Graph *graph, int radix)
{   
    printf("count sort for radix %d\n", radix);
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
    // auxillary arrays, allocated at the start up of the program
    int *vertex_count = (int *)malloc(graph->num_vertices * sizeof(int)); //needed for coutning sort

    #pragma omp parallel
    {
        int key;
        int pos;
        int P = omp_get_num_threads();
        int tid = omp_get_thread_num(); // used for parallization
        int offset_start, offset_end,vert_start, vert_end; // used for parallization
        int nbit = 8; // used for radix sort
        
        offset_start = tid * (graph->num_edges/P);
        vert_start = tid * (graph->num_vertices/P);
        if(tid == (P-1)){
            offset_end = offset_start + (graph->num_edges/P) + (graph->num_edges%P); // edge case for last group of indices
            vert_end = vert_start + (graph->num_vertices/P) + (graph->num_vertices%P);
        }else{
            offset_end = offset_start + (graph->num_edges/P);
            vert_end = vert_start + (graph->num_vertices/P);
        }
        if(tid == 0){
            printf("num_vertices:%d\n", graph->num_vertices);
            printf("num_edges:%d\n", graph->num_edges);
        }
        //printf("edges %d to %d\n", offset_start, offset_end);
        //printf("vertices %d to %d\n", vert_start, vert_end);

        for(int i = vert_start; i < vert_end; i++){
            vertex_count[i] = 0;
        }
        #pragma omp barrier
        // count occurrence of key: id of a source vertex
        printf("Started count occurence\n");
        for(int i = offset_start; i < offset_end; ++i)
        {
            key = graph->sorted_edges_array[i].src;
            int temp = ((key>>(radix * nbit)) & 0xffff); 
            vertex_count[temp]++;
            //printf("tid:%d vertex_count[%d]=%d\n", tid, i, vertex_count[i]);
            #pragma omp barrier
        }
        // transform to cumulative sum
        if(tid == 0){
            printf("tid:%d vertex_count[%d]=%d\n",tid, 0,vertex_count[0]);
            for(int i = 1; i < graph->num_vertices; i++)  {
                //printf("tid:%d vertex_count[%d]=%d\n",tid, i,vertex_count[i]);
                vertex_count[i] += vertex_count[i - 1];
            }
        }
        #pragma omp barrier
        // fill-in the sorted array of edges
        printf("Started fill-in\n");
        for(int i = offset_end - 1; i >= offset_start; --i)
        {
            key = graph->sorted_edges_array[i].src;
            int temp = (key >> (radix * nbit)) & 0xffff; // find 8-bit radix
            pos = vertex_count[temp]--;
            sorted_edges_array[pos] = graph->sorted_edges_array[i];
            //vertex_count[temp]--;
        }
        printf("tid(%d) finished\n",tid);
    } //end parallel region

    free(vertex_count);
    free(graph->sorted_edges_array);
    graph->sorted_edges_array = sorted_edges_array;
    printEdgeArray(graph->sorted_edges_array, graph->num_edges);
    return graph; 
}

#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{
    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
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
