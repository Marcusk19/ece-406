# Architecture Of Parallel Computers ECE 506 Projects and Material

This repository holds all machine problems, projects, and supplementary material for the ECE 506 course.

## Course Description

The need for parallel and massively parallel computers. Taxonomy of parallel computer architecture, and programming models for parallel architectures. Example parallel algorithms. Shared-memory vs. distributed-memory architectures. Correctness and performance issues. Cache coherence and memory consistency. Bus-based and scalable directory-based multiprocessors. Interconnection-network topologies and switch design. Brief overview of advanced topics such as multiprocessor pre-fetching and speculative parallel execution. Credit is not allowed for more than one course in this set: ECE 406, ECE 506.

## Please note the following

* **It is recommended to always to go back to the book as a reference.**
    * (NEW) Fundamentals of Parallel Multi-core Architecture (Chapman & Hall/CRC Computational Science) 1st Edition-2015 by Yan Solihin (Author).
    * Fundamentals of Parallel Computer Architecture Paperback-2009 by Yan Solihin (Author).
* For now this repository is for FALL 2021 Class any future updates are the responsibility of future TAs if they want to use it.
* This framework developed on Ubuntu-Linux. It is your responsibility to make it work on windows or other systems, though it will run on any Linux machine that has a terminal.

### Initial setup
1. Clone the repository. enter the directory and initialize the submodules
  ```bash
  git clone https://github.ncsu.edu/atmughra/ECE-506-406-Projects.git
  UnityID:
  PASSWord:
  cd ECE-506-406-Projects
  ```

## Project 1
### Benchmark format
* Edge lists are provided in textual format and they are formatted as source → destination separated by a tab.
*Src*  →  *Dest*

### Initial setup
1. After cloning the repository. enter the directory
  ```bash
  cd 01_Project-1/code
  ```
2. compile and run the code depending on which rule (openmp/mpi/hybrid) will generate a separate binary

3. run the code in openmp (bin/run-graph-openmp)
  ```bash
  make run-openmp
  ```
4. run the code in mpi (bin/run-graph-mpi)
  ```bash
  make run-mpi
  ```
5. run the code in Hyrbid (openmp+mpi) (bin/run-graph-hybrid)
  ```bash
  make run-hybrid
  ```

### Add Modifications
1. main.c/sort.c/sort.h files contain defines that are enabled depending on which rule you choose with the make file
  ```code
  make run-openmp

#ifdef OPENMP_HARNESS

// add code that you want to compile only for OpenMP

#endif

  ```
  ```code
  make run-mpi

#ifdef MPI_HARNESS

// add code that you want to compile only for MPI

#endif

  ```
  ```code
  make run-hybrid

#ifdef HYBRID_HARNESS

// add code that you want to compile only for MPI+OpenMP

#endif

  ```

## Organization
* `00_Architecture_of_Parallel_Computers`
  * `00_ARC_HPC_Guide` - How to setup ARC/HPC NCSU for future Projects/MachineProblems(MPs).
  * `01_OpenMP` - OpenMP tutorials and guides with examples within the slides.
  * `02_Slides` - Slides from chapter 1-12 from the old book 2009.
  * `03_papers` - Extra reading material and papers (optional).
  * `04_Fall_2018` - Updated slides and notes for the semester the class being taught.
* `Project-1` - Project 1
  * `code` - base code for the problem
    * `bin`- binaries are generated in this folder
    * `datasets` - sample data sets for the code
    * `include` - [.h] header files to be included in the code
        * `bfs.h` - push/pull breadth first search implementation for testing CSR correctness
        * `graph.h` - graph abstraction contains edgelists and inverse edge-lists with their mappings 
        * `bitmap.h` - bitmap implementation not multi-thread safe
        * `arrayQueue.h` - simple array based queue 
        * `edgelist.h`
        * `sort.h`  - count/radix sort implementation to be parallelized
        * `timer.h`
        * `vertex.h`
    * `obj` - [.o] object files generated here for linking
    * `src` - source code resides in this directory
        * `bfs.c` - push/pull breadth first search implementation for testing CSR correctness
        * `graph.c` - graph abstraction contains edge-lists and inverse edge-lists with their mappings
        * `bitmap.c` - bitmap implementation not multi-thread safe
        * `arrayQueue.c` - simple array based queue not multi-thread safe 
        * `edgelist.c`
        * `sort.c` - count/radix sort implementation to be parallelized
        * `timer.c`
        * `vertex.c`
        * `main.c` - **parallel/serial count/radix sort  called here**
    * `Makefile` - a file that compiles the base code
  * `doc` - Documents describes the Machine Problem

