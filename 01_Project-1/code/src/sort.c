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

struct Edge*countSortEdgesBySource (struct Edge *edgelist, int radix, int num_vertices, int num_edges, int *global_vertex_count)
{   
    // printf("count sort for radix %d\n", radix);
    // printf("num_edges=%d, num_vertices=%d\n", num_edges, num_vertices);
    int i;
    int key;
    int pos;
    int nbit = 8;
    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    struct Edge *sorted_edges_array = newEdgeArray(num_edges);

    int *vertex_count = (int*)malloc(num_vertices*sizeof(int)); // needed for counting sort
    // init vertex_count
    for(i = 0; i < num_vertices; i++){
        vertex_count[i] = 0;
    }
    // count occurence of key:
    for(i = 0; i < num_edges; i++){
        key = edgelist[i].src;
        int temp = (key>>(radix*nbit)) & 0xff;
        // printf("radix:%d, key:%d, temp:%d\n", radix, key, temp);
        vertex_count[temp]++;
        // printf("vertex_count[%d]=%d\n", temp, vertex_count[temp]);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce( vertex_count , global_vertex_count , num_vertices , MPI_INT , MPI_SUM , MPI_COMM_WORLD);
    // printf("radix(%d) count ok :)\n", radix);
    // transform to cumulative sum
    for(i = 1; i < num_vertices; i++){
        vertex_count[i] += vertex_count[i-1];
        //printf("vertex_count[%d]=%d\n", i, vertex_count[i]);
    }
    // printf("radix(%d) cumulative ok :)\n", radix);
    // printf("vertex_count[%d]=%d\n", num_vertices-1, vertex_count[num_vertices-1]);
    // fill-in sorted array
    for(i = num_edges-1; i >= 0; i--){
        key = edgelist[i].src;
        int temp  = (key>>(radix*nbit)) & 0xff;
        pos = vertex_count[temp] - 1;
        sorted_edges_array[pos] = edgelist[i];
        vertex_count[temp]--;
    }
    // printf("radix(%d) fill-in ok :)\n", radix);
    MPI_Barrier(MPI_COMM_WORLD);
    free(vertex_count);
    return sorted_edges_array; 
}




struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph)
{
    printf("*** Start Radix Sort Edges By Source MPI ***\n");
    #define nbit 8
    int myid, numprocs, offset_start, offset_end, edges_per_proc, num_vertices;
    int n = sizeof(graph->num_vertices);
    struct Edge *edge_buffer, *feedback_array;
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    int *global_vertex_count = (int*)malloc(graph->num_vertices*sizeof(int));

    for(int i = 0; i < graph->num_vertices; i++) global_vertex_count[i]=0;

    if(myid == 0) printf("n: %d\n",n);
    printf("I am %d of %d\n", myid,numprocs);
    
    // calculate how many edges per process
    if(myid == numprocs - 1){
        edges_per_proc = (graph->num_edges/numprocs) + (graph->num_edges%numprocs);
    } else {
        edges_per_proc = (graph->num_edges/numprocs);
    }
    // printf("edges_per_proc: %d\n", edges_per_proc);

    // create buffer to hold subset of edge array and feedback array
    edge_buffer = newEdgeArray(edges_per_proc);
    feedback_array = newEdgeArray(graph->num_edges);

    for(int i = 0; i < graph->num_edges; i++){
        feedback_array[i] = graph->sorted_edges_array[i];
    }

    num_vertices = graph->num_vertices;

    // initialize arrays for scatterv
    int scatter_displacements[numprocs];
    for(int i = 0; i < numprocs; i++){
        scatter_displacements[i] = edges_per_proc * i;
        // printf("rank(%d)scatter_displacements[%d]=%d\n",myid,i,scatter_displacements[i]);
    }
    int counts_send[numprocs];
    for(int i = 0; i < numprocs-1; i++){
        counts_send[i] = graph->num_edges/numprocs;
        // printf("counts_send[%d]=%d\n", i, counts_send[i]);
    }
    counts_send[numprocs-1] = (graph->num_edges/numprocs) + (graph->num_edges%numprocs);
    // printf("counts_send[%d]=%d\n", numprocs-1, counts_send[numprocs-1]);

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

    for(int radix = 0; radix < n; radix++){
        // scatter the edge array to each process
        if(myid == 0){
            MPI_Scatterv(feedback_array,  counts_send, scatter_displacements, edge_type , edge_buffer, edges_per_proc, edge_type , 0 , MPI_COMM_WORLD);
        } else {
            MPI_Scatterv( NULL , NULL , NULL , edge_type , edge_buffer , counts_send[myid] , edge_type , 0 , MPI_COMM_WORLD);
        }
        // process each local array
        edge_buffer = countSortEdgesBySource(edge_buffer, radix, num_vertices, edges_per_proc, global_vertex_count);
        MPI_Barrier( MPI_COMM_WORLD);
        // redistribute data among arrays
        if(numprocs == 1) continue; // rest of the loop is not needed for numprocs == 1
        struct Edge *window_buffer = newEdgeArray(graph->num_edges);
        MPI_Win window;
        MPI_Win_create(window_buffer, graph->num_edges*sizeof(window_buffer[0]), sizeof(window_buffer[0]), MPI_INFO_NULL, MPI_COMM_WORLD, &window);
        MPI_Win_fence(0, window);
        // printf("rank(%d) window created\n", myid);
        MPI_Barrier(MPI_COMM_WORLD);

        for(int i = 0; i < edges_per_proc; i++){
            int offset = (myid * (graph->num_edges/numprocs)) + i;
            if(myid != 0){
                struct Edge to_send = edge_buffer[i];
                MPI_Put(&to_send, 1, edge_type, 0, offset, 1, edge_type, window);
            } else {
                window_buffer[offset] = edge_buffer[i];
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        

        if(myid == 0){
            struct Edge *tmp_edges = newEdgeArray(graph->num_edges);
            int tmp_index = 0;
            for(int vertex = 0; vertex < num_vertices; vertex++){
                int current_radix = (vertex >> (radix * nbit)) & 0xff;
                if(global_vertex_count[current_radix] == 0) continue;
                int num_keys = global_vertex_count[current_radix];

                for(int rank = 0; rank < numprocs; rank++){
                    if(num_keys <= 0) break;
                    int sub_len;
                    if(rank == numprocs - 1){
                        sub_len = edges_per_proc + (graph->num_edges%numprocs);
                    } else {
                        sub_len = edges_per_proc;
                    }

                    for(int i = 0; i < sub_len; i++){
                        int local_index = (rank * (graph->num_edges/numprocs)) + i;
                        int key = window_buffer[local_index].src;
                        int temp = (key >> (radix*nbit)) & 0xff;
                        if(temp == current_radix && num_keys > 0){
                            // printf("matched temp %d to current_radix %d\n", temp, current_radix);
                            // printf("radix(%d), rank(%d), tmp_index = %d\n",radix, rank, tmp_index);
                            tmp_edges[tmp_index] = window_buffer[local_index];
                            num_keys--;
                            tmp_index++;
                        }
                        if(num_keys <= 0) break;
                    }
                    if(num_keys <= 0) break;
                }   
                if(tmp_index >= graph->num_edges) break;
            }
            for(int i = 0; i < graph->num_edges; i++){
                feedback_array[i] = tmp_edges[i];
            }
            printf("radix(%d) swapped edgelist\n", radix);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_free(&window); 
    } // end of for radix
    
    if(myid == 0){
        // swap arrays
        struct Edge *temp_pointer = graph->sorted_edges_array;
        graph->sorted_edges_array = feedback_array;
        feedback_array = temp_pointer;
    }
    free(global_vertex_count);
    printf("rank(%d) freed global_vertex_count\n", myid);
    // MPI_Bcast( graph->sorted_edges_array , graph->num_edges , edge_type , 0 , MPI_COMM_WORLD); // broadcast sorted graph to all processors

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
