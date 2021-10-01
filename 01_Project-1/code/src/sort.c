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

#ifdef OPENMP_HARNESS
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
        //printf("Started count occurence\n");
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
            pos = vertex_count[tid*graph->num_vertices+temp]-1; // here a specific thread's index = tid*graph->num_vertices+temp
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

struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph)
{
    printf("*** START Radix Sort Edges By Source OpenMP *** \n");
    int n = sizeof(graph->num_vertices);
    printf("n:%d \n", n);

    
    for(int i = 0; i < n; i++){
        graph = countSortEdgesBySource(graph, i);
    }

    
    printEdgeArray(graph->sorted_edges_array, graph->num_edges);
    return graph;
}
#endif

//-------------------- MPI_HARNESS----------------------------------
#ifdef MPI_HARNESS
struct Edge*countSortEdgesBySource (struct Edge *edgelist, int radix, int num_vertices, int num_edges)
{   
    printf("count sort for radix %d\n", radix);
    printf("num_edges=%d, num_vertices=%d\n", num_edges, num_vertices);
    int i;
    int key;
    int pos;
    int nbit = 8;
    struct Edge *sorted_edges_array = newEdgeArray(num_edges);
    int *vertex_count = (int*)malloc(num_vertices*sizeof(int)); // needed for counting sort
    // init vertex_count
    for(i = 0; i < num_vertices; i++){
        vertex_count[i] = 0;
    }
    printf("radix(%d) init. ok :) \n", radix);
    // count occurence of key:
    for(i = 0; i < num_edges; i++){
        key = edgelist[i].src;
        int temp = (key>>(radix*nbit)) & 0xff;
        //printf("radix:%d, key:%d, temp:%d\n", radix, key, temp);
        vertex_count[temp]++;
        //printf("vertex_count[%d]=%d\n", temp, vertex_count[temp]);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printf("radix(%d) count ok :)\n", radix);
    // transform to cumulative sum
    for(i = 1; i < num_vertices; i++){
        vertex_count[i] += vertex_count[i-1];
        //printf("vertex_count[%d]=%d\n", i, vertex_count[i]);
    }
    printf("radix(%d) cumulative ok :)\n", radix);
    // printf("vertex_count[%d]=%d\n", num_vertices-1, vertex_count[num_vertices-1]);
    // fill-in sorted array
    for(i = num_edges-1; i >= 0; i--){
        key = edgelist[i].src;
        int temp  = (key>>(radix*nbit)) & 0xff;
        pos = vertex_count[temp] - 1;
        sorted_edges_array[pos] = edgelist[i];
        vertex_count[temp]--;
    }
    printf("radix(%d) fill-in ok :)\n", radix);
    MPI_Barrier(MPI_COMM_WORLD);
    free(vertex_count);
    // printEdgeArray(sorted_edges_array, num_edges);
    return sorted_edges_array; 
}

struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph)
{
    printf("*** Start Radix Sort Edges By Source MPI ***\n");
    int myid, numprocs, offset_start, offset_end, edges_per_proc, num_vertices;
    int n = sizeof(graph->num_vertices);

    struct Edge *my_edgelist;
    struct Edge *tmp_edges;
    struct Edge *edge_buffer;
    

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    if(myid == 0) printf("n: %d\n",n);
    printf("I am %d of %d\n", myid,numprocs);
    

    edges_per_proc = (graph->num_edges/numprocs);
    printf("edges_per_proc: %d\n", edges_per_proc);
    // create buffer to hold subset of edge array
    if(myid == numprocs - 1){
        edges_per_proc += graph->num_edges%numprocs;
    }
    edge_buffer = newEdgeArray(edges_per_proc);

    num_vertices = (graph->num_vertices/numprocs);

    int scatter_displacements[numprocs];
    for(int i = 0; i < numprocs; i++){
        scatter_displacements[i] = edges_per_proc * i;
        printf("rank(%d)scatter_displacements[%d]=%d\n",myid,i,scatter_displacements[i]);
    }
    int counts_send[numprocs];
    for(int i = 0; i < numprocs-1; i++){
        counts_send[i] = edges_per_proc;
        printf("counts_send[%d]=%d\n", i, counts_send[i]);
    }
    counts_send[numprocs-1] = (graph->num_edges/numprocs) + (graph->num_edges%numprocs);
    printf("counts_send[%d]=%d\n", numprocs-1, counts_send[numprocs-1]);

    // build datatype for struct Edge
        int lengths[2] = {1, 1};
        MPI_Datatype edge_type;
        MPI_Aint displacements[2];
        struct Edge dummy_edge;
        MPI_Aint base_address;
        MPI_Get_address(&dummy_edge, &base_address);
        MPI_Get_address( &dummy_edge.src , &displacements[0]);
        MPI_Get_address( &dummy_edge.dest , &displacements[1]);
        displacements[0] = MPI_Aint_diff(displacements[0], base_address);
        displacements[1] = MPI_Aint_diff(displacements[1], base_address);
        MPI_Datatype types[2] = {MPI_INT, MPI_INT};
        MPI_Type_create_struct(2, lengths, displacements, types, &edge_type);
        MPI_Type_commit(&edge_type);
    // split data to processors using MPI_Scatter
    
    for(int radix = 0; radix < 1; radix++){
        if(myid == 0){
            MPI_Scatterv(graph->sorted_edges_array,  counts_send, scatter_displacements, edge_type , edge_buffer, edges_per_proc, edge_type , 0 , MPI_COMM_WORLD);
        } else {
            MPI_Scatterv( NULL , NULL , NULL , edge_type , edge_buffer , counts_send[myid] , edge_type , 0 , MPI_COMM_WORLD);
        }
        printf("rank(%d) received from scatterv\n", myid);
        for(int i = edges_per_proc - 1; i >= edges_per_proc - 10; i--){
            printf("rank(%d)edge_buffer[%d]=%d->%d\n", myid, i, edge_buffer[i].src, edge_buffer[i].dest);
        }
        edge_buffer = countSortEdgesBySource(edge_buffer, radix, num_vertices, edges_per_proc);
        printf("rank(%d) count sorted radix(%d)\n", myid, radix);

        // if(myid == 1){
        //     printEdgeArray(edge_buffer, edges_per_proc);
        // }


        if(myid == 0){
            my_edgelist = newEdgeArray(graph->num_edges);
            MPI_Gatherv(edge_buffer , edges_per_proc , edge_type , my_edgelist , counts_send, scatter_displacements, edge_type , 0 , MPI_COMM_WORLD);
            free(graph->sorted_edges_array);
            graph->sorted_edges_array = my_edgelist;
            //if(radix == 0)printEdgeArray(graph->sorted_edges_array, graph->num_edges);
        //     free(my_edgelist);
        } else {
            MPI_Gatherv(edge_buffer, edges_per_proc, edge_type, NULL, NULL, NULL, edge_type, 0, MPI_COMM_WORLD);
        }

        // merge arrays
        
    }


    for(int i = graph->num_edges-1; i >= graph->num_edges-11; i--){
        printf("%d->%d\n", my_edgelist[i].src, my_edgelist[i].dest);
    }
    

    
    // MPI_Barrier(MPI_COMM_WORLD);
    // printf("--------rank(%d) edgelist--------\n", myid);
    // printEdgeArray(my_edgelist, 10);
    
    // for(int i = 0; i < n; i++){
    //     my_edgelist = countSortEdgesBySource(graph->sorted_edges_array, i, graph->num_vertices, graph->num_edges);
    // }
    //free(graph->sorted_edges_array);
    //graph->sorted_edges_array = my_edgelist;
    // printEdgeArray(graph->sorted_edges_array, graph->num_edges);
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
