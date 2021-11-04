/*******************************************************
                          main.cc
********************************************************/

#include <stdlib.h>
#include <assert.h>
#include <fstream>
using namespace std;

#include "cache.h"

int main(int argc, char *argv[])
{
	printf("===== 506 Personal information =====\n");
	printf("Marcus Kok\n");
	printf("mskok\n");
	printf("undergrad section? YES\n");
	ifstream fin;
	FILE * pFile;

	if(argv[1] == NULL){
		 printf("input format: ");
		 printf("./smp_cache <cache_size> <assoc> <block_size> <num_processors> <protocol> <trace_file> \n");
		 exit(0);
        }

	/*****uncomment the next five lines*****/
	int cache_size = atoi(argv[1]);
	int cache_assoc= atoi(argv[2]);
	int blk_size   = atoi(argv[3]);
	int num_processors = atoi(argv[4]);/*1, 2, 4, 8*/
	int protocol   = atoi(argv[5]);	 /*0:MSI, 1:MESI, 2:Dragon*/

	char *fname =  (char *)malloc(20);
 	fname = argv[6];
	char line[256];
	
	//****************************************************//
	printf("===== 506 SMP Simulator configuration =====\n");
	//*******print out simulator configuration here*******//
	//****************************************************//
	printf("L1_SIZE: %d\n", cache_size);
	printf("L1_ASSOC: %d\n", cache_assoc);
	printf("L1_BLOCKSIZE: %d\n", blk_size);
	printf("NUMBER OF PROCESSORS: %d\n", num_processors);
	if(protocol == MSI) printf("COHERENCE PROTOCOL: MSI\n");
	else printf("COHERENCE PROTOCOL: MESI\n");
	printf("TRACE FILE: %s\n", fname);
 
	//*********************************************//
	//*****create an array of caches here**********//
	//*********************************************//	
	Cache* cachesArray[num_processors];
	for(int i = 0; i < num_processors; i++){
		cachesArray[i] = new Cache(cache_size, cache_assoc, blk_size);
	}

	pFile = fopen (fname,"r");
	if(pFile == 0)
	{   
		printf("Trace file problem\n");
		exit(0);
	}
	///******************************************************************//
	//**read trace file,line by line,each(processor#,operation,address)**//
	//*****propagate each request down through memory hierarchy**********//
	//*****by calling cachesArray[processor#]->Access(...)***************//
	///******************************************************************//
	while(fgets(line, sizeof(line), pFile)){
		string currentLine(line); // making copy of char* line for substr()
		// parse out the processor, operation, and address from line
		int proc = line[0] - 48; // ascii conversion
		uchar op = (uchar)line[2];
		string parsedAddr = currentLine.substr(4, currentLine.length());
		ulong addr = strtoul(parsedAddr.c_str(), nullptr, 16); // convert string into ulong value for addr
		cachesArray[proc]->Access(addr, op, protocol, cachesArray, num_processors, proc); // call access function
		// snoop all other caches
		for(int i = 0; i < num_processors; i++){
			if(i != proc) cachesArray[i]->Snoop(addr, op, protocol);
		}


	}
	fclose(pFile);

	//********************************//
	//print out all caches' statistics //
	//********************************/
	
	for(int i = 0; i < num_processors; i++){
		printf("============ Simulation results (Cache %d) ============\n", i);
		/****print out the rest of statistics here.****/
		cachesArray[i]->printStats();
	}
}
