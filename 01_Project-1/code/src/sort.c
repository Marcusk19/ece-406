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

struct Graph *countSortEdgesBySource (struct Graph *graph)
{


    int printThreads = 0;
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);
    #pragma omp parallel
    {
        int i;
        int key;
        int pos;
        int P;
        int tid = omp_get_thread_num();
        int offset_start, offset_end;

        P = omp_get_num_threads();
        int *vertex_count = (int*)malloc(P*graph->num_edges*sizeof(int));
        printf("malloc:%d\n", graph->num_edges*sizeof(int));
        
        #pragma omp barrier

        offset_start = tid * (graph->num_edges/P);

        if(tid == (P-1)){
            offset_end = offset_start + (graph->num_edges/P) + (graph->num_edges%P);
        }else{
            offset_end = offset_start + (graph->num_edges/P);
        }
        printf("num_edges:%d\n", graph->num_edges);
        printf("%d to %d\n", offset_start, offset_end);
        #pragma omp barrier

        // auxiliary arrays, allocated at the start up of the program
        //int *vertex_count = (int *)malloc(graph->num_vertices * sizeof(int)); // needed for Counting Sort

        printf("Started array initialization\n");
        for(i = offset_start; i < offset_end; ++i)
        {
            //printf("%d\n", i);
            vertex_count[i] = 0;
            if(printThreads)printf("Initialize array iteration=%d thread=%d\n",i,(int)pthread_self());
        }
        printf("Finished array initialization\n");
        #pragma omp barrier

        // count occurrence of key: id of a source vertex
        printf("Started count occurence\n");
        for(i = offset_start; i < offset_end; ++i)
        {
            key = graph->sorted_edges_array[i].src;
            //if(printThreads)printf("Count occurence iteration=%d thread=%d\n",i,(int)pthread_self());
            vertex_count[key]++;
            if(printThreads)printf("Count occurence iteration=%d thread=%d\n",i,(int)pthread_self());
        }
        printf("Finished count occurence\n");
        #pragma omp barrier

        // transform to cumulative sum
        printf("Started cumulative sum\n");
        for(i = 1; i < graph->num_vertices; ++i)
        {
            vertex_count[i] += vertex_count[i - 1];
            if(printThreads)printf("Cumulative sum iteration=%d thread=%d\n",i,(int)pthread_self());
        }
        printf("Finished cumulative sum\n");
        #pragma omp barrier

        // fill-in the sorted array of edges
        printf("Started fill-in\n");
        for(i = offset_end - 1; i >= offset_start; --i)
        {
            key = graph->sorted_edges_array[i].src;
            pos = vertex_count[key] - 1;
            sorted_edges_array[pos] = graph->sorted_edges_array[i];
            vertex_count[key]--;
            if(printThreads)printf("Fill in sorted array iteration=%d thread=%d\n",i,(int)pthread_self());
        }
        printf("Finished fill-in\n");
        #pragma omp barrier
        free(vertex_count);
    }  
    free(graph->sorted_edges_array);
    graph->sorted_edges_array = sorted_edges_array;
    printEdgeArray(sorted_edges_array, graph->num_edges);
    return graph; 

}

#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{

    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
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