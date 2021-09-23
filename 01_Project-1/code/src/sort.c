#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
#include <pthread.h>
// #include <mpi.h>

#include "sort.h"
#include "graph.h"

struct Graph* countSortEdgesBySource (struct Graph* graph){
    //printf("--Parallized Count Sort-- \n");
    int printThreads = 0;
    //int i = 0;
    int key;
    int pos;
    int P;
    struct Edge *sorted_edges_array = newEdgeArray(graph->num_edges);

    // auxiliary arrays, allocated at the start up of the program
    //int *vertex_count = (int*)malloc(graph->num_vertices*sizeof(int)); // needed for Counting Sort
    
    // from piazza
    printf("num_edges = %d\n", graph->num_edges);

    printf("parallel region start \n");
    #pragma omp parallel
    {
        // start parallel region
        int i = 0;
        int tid = omp_get_thread_num();
        printf("tid = %d\n", tid);


            P = omp_get_num_threads();
            //printf("P = %d\n", P);
            int *vertex_count = (int*)malloc(P*(graph->num_vertices * sizeof(int)));
        

        #pragma omp barrier
        int offset_start = tid * (graph->num_edges / P);
        int offset_end;
        if(tid == (P-1)){
            offset_end = offset_start + (graph->num_edges/P) + (graph->num_edges%P);
        }
        else{
            offset_end = offset_start + (graph->num_edges/P);
        }

        printf("offset_start = %d and offset_end = %d\n", offset_start, offset_end);
        // this loop will excecute in parallel but need to split i = offset_start, i < offset_end
        // according to tid
        #pragma omp barrier

        printf("Start initialization\n");
        for(int i = offset_start; i < offset_end; i++){
            vertex_count[i] = 0; 
            if(printThreads)printf("Initialize array iteration=%d thread=%d\n", i, (int)pthread_self());
        }
        #pragma omp barrier

        printf("\nStart count\n");
        for(int i = offset_start; i < offset_end; i++){
            key = graph->sorted_edges_array[i].src;
            vertex_count[key]++;
            if(printThreads)printf("Count occurence iteration=%d thread=%d\n",i,(int)pthread_self());
        }
        #pragma omp barrier

        printf("Start cumulative sum\n");
        for(i = offset_start; i < graph->num_vertices; i++){
            vertex_count[i] += vertex_count[i - 1];
            if(printThreads)printf("Cumulative sum iteration=%d thread=%d\n",i,(int)pthread_self());
        }
        #pragma omp barrier

        printf("Start array fill-in\n");
        for(i = offset_start; i < offset_end; i++){
            key = graph->sorted_edges_array[i].src;
            pos = vertex_count[key] - 1;
            sorted_edges_array[pos] = graph->sorted_edges_array[i];
            vertex_count[key]--;
            if(printThreads)printf("Fill in sorted array iteration=%d thread=%d\n",i,(int)pthread_self());
        }
        free(vertex_count);

    } // end parallel region

 
    printf("Freed vertex_count\n");
    free(graph->sorted_edges_array);

    graph->sorted_edges_array = sorted_edges_array;
    //printEdgeArray(sorted_edges_array, graph->num_edges);
    return graph;
    
    
    /*
    #pragma omp parallel for shared(vertex_count,i)
    for(i=0; i < graph->num_vertices; ++i) {
        vertex_count[i] = 0;
        if(printThreads)printf("Initialize array iteration=%d thread=%d\n",i,(int)pthread_self());
    }
    #pragma omp barrier

    // count occurrence of key: id of a source vertex
    #pragma omp parallel for shared(vertex_count,key,sorted_edges_array,i)
    for(i = 0; i < graph->num_edges; ++i) {
        key = graph->sorted_edges_array[i].src;
        //if(printThreads)printf("Count occurence iteration=%d thread=%d\n",i,(int)pthread_self());
        vertex_count[key]++;
        if(printThreads)printf("Count occurence iteration=%d thread=%d\n",i,(int)pthread_self());
    }
    #pragma omp barrier

    #pragma omp parallel for shared(vertex_count,i)
    // transform to cumulative sum
    for(i = 1; i < graph->num_vertices; ++i) {
        vertex_count[i] += vertex_count[i - 1];
        if(printThreads)printf("Cumulative sum iteration=%d thread=%d\n",i,(int)pthread_self());
    }
    #pragma omp barrier

    #pragma omp parallel for shared(vertex_count,key,pos,sorted_edges_array,i)
    // fill-in the sorted array of edges
    for(i = graph->num_edges - 1; i >= 0; --i) {
        key = graph->sorted_edges_array[i].src;
        pos = vertex_count[key] - 1;
        sorted_edges_array[pos] = graph->sorted_edges_array[i];
        vertex_count[key]--;
        if(printThreads)printf("Fill in sorted array iteration=%d thread=%d\n",i,(int)pthread_self());
    }
    #pragma omp barrier

    free(vertex_count);
    free(graph->sorted_edges_array);

    graph->sorted_edges_array = sorted_edges_array;
    printEdgeArray(sorted_edges_array, graph->num_edges);
    return graph; 
    */

}

struct Graph* radixSortEdgesBySourceOpenMP (struct Graph* graph){

    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
    
    return NULL;
} 
// struct Graph* radixSortEdgesBySourceMPI (struct Graph* graph){

//         // printf("*** START Radix Sort Edges By Source MPI*** \n");

// }
// struct Graph* radixSortEdgesBySourceHybrid (struct Graph* graph){

//         // printf("*** START Radix Sort Edges By Source Hybrid*** \n");

// }
