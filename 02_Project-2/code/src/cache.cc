/*******************************************************
                          cache.cc
********************************************************/

#include <stdlib.h>
#include <assert.h>
#include "cache.h"
using namespace std;

Cache::Cache(int s,int a,int b )
{
   ulong i, j;
   reads = readMisses = writes = 0; 
   writeMisses = writeBacks = currentCycle = 0;

   size       = (ulong)(s);
   lineSize   = (ulong)(b);
   assoc      = (ulong)(a);   
   sets       = (ulong)((s/b)/a);
   numLines   = (ulong)(s/b);
   log2Sets   = (ulong)(log2(sets));   
   log2Blk    = (ulong)(log2(b));   
  
   //*******************//
   //initialize your counters here//
   //*******************//
   busRdX = invalidations = c2cTransfers = flushes = interventions = memoryTransactions = 0;
 
   tagMask =0;
   for(i=0;i<log2Sets;i++)
   {
      tagMask <<= 1;
        tagMask |= 1;
   }
   /**create a two dimentional cache, sized as cache[sets][assoc]**/ 
   cache = new cacheLine*[sets];
   for(i=0; i<sets; i++)
   {
      cache[i] = new cacheLine[assoc];
      for(j=0; j<assoc; j++) 
      {
      cache[i][j].invalidate();
      }
   }      
}

/**you might add other parameters to Access()
since this function is an entry point 
to the memory hierarchy (i.e. caches)**/
void Cache::Access(ulong addr,uchar op, int protocol)
{
   currentCycle++;/*per cache global counter to maintain LRU order 
         among cache ways, updated on every cache access*/

   cacheLine *line;
   switch(protocol){
      case MSI:
         if(op == 'w') writes++;
         else reads++;

         line = findLine(addr);
         if(line == NULL){
            // miss
            memoryTransactions++;
            if(op == 'w'){
               writeMisses++;
            }
            else readMisses++;

            cacheLine *newline = fillLine(addr);
            if(op == 'w'){
               newline->setFlags(DIRTY); // I->M
               busRdX++; // read with intention to write
            }
            else newline->setFlags(VALID); // I->S
         }
         else{
            // hit, update LRU and update dirty flag
            updateLRU(line);
            if(line->getFlags() == VALID){
               if(op == 'w'){
                  line->setFlags(DIRTY); // S->M
                  memoryTransactions++;
                  busRdX++;
               }
               else line->setFlags(VALID); // S->S
            }
            else if(line->getFlags() == DIRTY){
               line->setFlags(DIRTY); // M->M
            }
         }
      break;
      //================== MESI =======================
      case MESI:
         if(op == 'w') writes++;
         else reads++;

         line = findLine(addr);
         if(line == NULL){
            // miss
            memoryTransactions++;
            if(op == 'w') writeMisses++;
            else readMisses++;

            cacheLine *newLine = fillLine(addr);
            if(op == 'w'){
               newLine->setFlags(DIRTY); // I->M
               busRdX++;
            }
            else if(op == 'r' && newLine->copyExists()) newLine->setFlags(VALID); // I->S
            else if(op == 'r' && !(newLine->copyExists())) newLine->setFlags(EXCLUSIVE); // I->E

         }
         else {
            // hit, update LRU
            updateLRU(line);
            if(line->getFlags() == EXCLUSIVE){
               if(op == 'w'){
                  line->setFlags(DIRTY); // E->M
                  line->setCopyFlag(true);
               }
            }
            if(line->getFlags() == VALID){
               if(op == 'w'){
                  line->setFlags(DIRTY); // S->M
                  line->setCopyFlag(true);
               }
            }
            // by default M->M
         }
      break;
   }
}

void Cache::Snoop(ulong addr, uchar op, int protocol){

   cacheLine* line;

   switch (protocol){
      case MSI:
         line = findLine(addr);
         if(line != NULL){
            // hit (not invalid state)
            // updateLRU(line);

            if(line->getFlags() == VALID){
               if(op == 'w'){
                  line->invalidate(); // S->I
                  invalidations++;
               }
               // S->S otherwise
            }
            else if(line->getFlags() == DIRTY){
               writeBacks++; // dirty flag so need to do a writeback
               memoryTransactions++; // memory transaction from writeback
               if(op == 'w'){
                  line->invalidate(); // M->I
                  invalidations++;
                  flushes++;
               }
               else {
                  line->setFlags(VALID); // M->S
                  flushes++;
                  interventions++; // intervention from moving to shared state
               }
            }
         }
      break;

      case MESI:
         line = findLine(addr);
         if(line !=  NULL){
            // line exists in cache
            if(line->getFlags() == VALID){
               if(op == 'w'){
                  line->invalidate(); // S->I
                  invalidations++;
               } // otherwise S->S
            }
            else if(line->getFlags() == EXCLUSIVE){
               if(op == 'w'){
                  line->invalidate(); // E->I
                  invalidations++;
               }
               else{
                  line->setFlags(VALID); // E->S
                  line->setCopyFlag(true);
                  interventions++;
               }
            }
            else if(line->getFlags() == DIRTY){
               if(op == 'w'){
                  line->invalidate(); // M->I
                  invalidations++;
                  flushes++;
               }
               else{
                  line->setFlags(VALID); // M->S
                  flushes++;
                  interventions++;
               }
            }
            // otherwise I->I
         }
      break;
   }
}

/*look up line*/
cacheLine * Cache::findLine(ulong addr)
{
   ulong i, j, tag, pos;
   
   pos = assoc;
   tag = calcTag(addr);
   i   = calcIndex(addr);
  
   for(j=0; j<assoc; j++)
   if(cache[i][j].isValid())
           if(cache[i][j].getTag() == tag)
      {
           pos = j; break; 
      }
   if(pos == assoc)
   return NULL;
   else
   return &(cache[i][pos]); 
}

/*upgrade LRU line to be MRU line*/
void Cache::updateLRU(cacheLine *line)
{
  line->setSeq(currentCycle);  
}

/*return an invalid line as LRU, if any, otherwise return LRU line*/
cacheLine * Cache::getLRU(ulong addr)
{
   ulong i, j, victim, min;

   victim = assoc;
   min    = currentCycle;
   i      = calcIndex(addr);
   
   for(j=0;j<assoc;j++)
   {
      if(cache[i][j].isValid() == 0) return &(cache[i][j]);     
   }  
   for(j=0;j<assoc;j++)
   {
    if(cache[i][j].getSeq() <= min) { victim = j; min = cache[i][j].getSeq();}
   } 
   assert(victim != assoc);
   
   return &(cache[i][victim]);
}

/*find a victim, move it to MRU position*/
cacheLine *Cache::findLineToReplace(ulong addr)
{
   cacheLine * victim = getLRU(addr);
   updateLRU(victim);
  
   return (victim);
}

/*allocate a new line*/
cacheLine *Cache::fillLine(ulong addr)
{ 
   ulong tag;
  
   cacheLine *victim = findLineToReplace(addr);
   assert(victim != 0);
   if(victim->getFlags() == DIRTY){
       writeBack(addr); // function increments writeBacks counter
       memoryTransactions++; // memory transaction as part of writeback
   }
       
   tag = calcTag(addr);   
   victim->setTag(tag);
   victim->setFlags(VALID);    
   /**note that this cache line has been already 
      upgraded to MRU in the previous function (findLineToReplace)**/

   return victim;
}

void Cache::printStats()
{ 
   double missRate = ((double)writeMisses + (double)readMisses) / ((double)reads + (double)writes) * 100;
   printf("01. number of reads:                    %ld\n", reads);
   printf("02. number of read misses:              %ld\n", readMisses);
   printf("03. number of writes:                   %ld\n", writes);
   printf("04. number of write misses:             %ld\n", writeMisses);
   printf("05. total miss rate:                    %.2f%%\n", missRate);
   printf("06. number of writebacks:               %ld\n", writeBacks);
   printf("07. number of cache-to-cache transfers: %ld\n", c2cTransfers);
   printf("08. number of memory transactions:      %ld\n", memoryTransactions);
   printf("09. number of interventions:            %ld\n", interventions);
   printf("10. number of invalidations:            %ld\n", invalidations);
   printf("11. number of flushes:                  %ld\n", flushes);
   printf("12. number of BusRdX:                   %ld\n", busRdX);
   /****follow the ouput file format**************/
}
