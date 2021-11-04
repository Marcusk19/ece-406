/*******************************************************
                          cache.h
********************************************************/

#ifndef CACHE_H
#define CACHE_H

#include <cmath>
#include <iostream>

typedef unsigned long ulong;
typedef unsigned char uchar;
typedef unsigned int uint;

#define MSI 0
#define MESI 1

/****add new states, based on the protocol****/
enum{
	INVALID = 0,
	VALID,
	DIRTY,
   // added exclusive for MESI
   EXCLUSIVE
};

class cacheLine 
{
protected:
   ulong tag;
   ulong Flags;   // 0:invalid, 1:valid, 2:dirty 
   ulong seq; 
   bool C; // flag for MESI
 
public:
   cacheLine()            { tag = 0; Flags = 0; C = false;}
   ulong getTag()         { return tag; }
   ulong getFlags()			{ return Flags;}
   ulong getSeq()         { return seq; }
   void setSeq(ulong Seq)			{ seq = Seq;}
   void setFlags(ulong flags)			{  Flags = flags;}
   void setTag(ulong a)   { tag = a; }
   void invalidate()      { tag = 0; Flags = INVALID; }//useful function
   bool isValid()         { return ((Flags) != INVALID); }

   bool copyExists() { return C;} // for MESI
   void setCopyFlag(bool set){ C = set; }
};

class Cache
{
protected:
   ulong size, lineSize, assoc, sets, log2Sets, log2Blk, tagMask, numLines;
   ulong reads,readMisses,writes,writeMisses,writeBacks;

   //******///

   //add coherence counters here///
   //******///
   ulong busRdX, invalidations, c2cTransfers, flushes, interventions, memoryTransactions;
   // ulong dirtyBlocks, c2cTransfers, memoryTransactions, interventions,
   //       invalidations, flushes, busRdX;

   cacheLine **cache;
   ulong calcTag(ulong addr)     { return (addr >> (log2Blk) );}
   ulong calcIndex(ulong addr)  { return ((addr >> log2Blk) & tagMask);}
   ulong calcAddr4Tag(ulong tag)   { return (tag << (log2Blk));}
   
public:
    ulong currentCycle;  
     
    Cache(int,int,int);
   ~Cache() { delete cache;}
   
   cacheLine *findLineToReplace(ulong addr);
   cacheLine *fillLine(ulong addr);
   cacheLine * findLine(ulong addr);
   cacheLine * getLRU(ulong);
   
   ulong getRM(){return readMisses;} ulong getWM(){return writeMisses;} 
   ulong getReads(){return reads;}ulong getWrites(){return writes;}
   ulong getWB(){return writeBacks;}
   
   void writeBack(ulong)   {writeBacks++;}
   void Access(ulong,uchar,int, Cache**, int, int);
   void printStats();
   void updateLRU(cacheLine *);

   //******///
   //add other functions to handle bus transactions///
   //******///
   void Snoop(ulong, uchar, int);
   bool isCopy(ulong, Cache**, int, int);

};

#endif
