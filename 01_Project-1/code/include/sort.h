#ifndef SORT_H
#define SORT_H

#include "edgelist.h"
#include "graph.h"
// Order edges by id of a source vertex,
// using the Counting Sort
// Complexity: O(E + V)

<<<<<<< HEAD
// struct Graph* radixSortEdgesBySourceOpenMP (struct Graph* graph);
// struct Graph* radixSortEdgesBySourceMPI (struct Graph* graph);
// struct Graph* radixSortEdgesBySourceHybrid (struct Graph* graph);
struct Graph* countSortEdgesBySource (struct Graph* graph, int nbit);
=======
#ifdef OPENMP_HARNESS
struct Graph *radixSortEdgesBySourceOpenMP (struct Graph *graph);
#endif

#ifdef MPI_HARNESS
struct Graph *radixSortEdgesBySourceMPI (struct Graph *graph);
#endif

#ifdef HYBRID_HARNESS
struct Graph *radixSortEdgesBySourceHybrid (struct Graph *graph);
#endif

struct Graph *countSortEdgesBySource (struct Graph *graph, int radix);
>>>>>>> revert

extern int numThreads;

#endif