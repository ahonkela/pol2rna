#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>

/*
#define mymalloc malloc
#define myrealloc realloc
#define myfree free
#define myprintf printf
*/


#include "mex.h"
#include "matrix.h"

#define mymalloc mxMalloc
#define myrealloc mxRealloc
#define myfree mxFree
#define myprintf mexPrintf



#define MAX_LINELENGTH 1024
#define n_linestoskip 4

#define n_chromosomes 25
const char *chromosomenames[]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT"};


#define index_readstart 0
#define index_readend 1
#define index_readstrand 2
#define index_readscore 3
#define index_nreads 4
#define index_maxreads 5
#define index_linenumber 6


void ***read_mappingfile(char *filename)
{
  FILE *f1;
  int mappings_available = 0;
  char tline1[MAX_LINELENGTH];
  char tempbuffer[MAX_LINELENGTH];
  char *readsuccess;

  int i, k, i1, i2;
  int chr_index;
  int nlines;

  void ***mappings;
  int nreads;
  int maxreads;
  int *readstarts;
  int *readends;
  signed char *readstrands;
  float *readscores;
  int *linenumbers;

  int line_chrindex;
  int line_readstart;
  int line_readend;
  signed char line_readstrand;
  float line_readscore;


  /*------------------------------
  Initialize data cell array
  ------------------------------*/
  myprintf("Initializing cell array\n",filename);

  mappings=(void ***)mymalloc(n_chromosomes*sizeof(void **));
  for (chr_index=0;chr_index<n_chromosomes;chr_index++)
  {
    mappings[chr_index] = (void **)mymalloc(7*sizeof(void *));

    mappings[chr_index][index_readstart] = NULL;
    mappings[chr_index][index_readend] = NULL;
    mappings[chr_index][index_readstrand] = NULL;
    mappings[chr_index][index_readscore] = NULL;
    mappings[chr_index][index_linenumber] = NULL;
    mappings[chr_index][index_nreads] = mymalloc(1*sizeof(int));
    mappings[chr_index][index_maxreads] = mymalloc(1*sizeof(int));
    *((int *)(mappings[chr_index][index_nreads])) = 0;
    *((int *)(mappings[chr_index][index_maxreads])) = 0;
  }


  /*------------------------------
    Check the file exists
    ------------------------------*/
  myprintf("Checking file existence\n",filename);
  f1 = fopen(filename,"rb");
  if (!f1) 
  {
    mappings_available=0; 
    myprintf("Could not open file [%s] for reading\n",filename);
    return(mappings);
  }


  /*------------------------------
    Skip first lines (header)
    ------------------------------*/
  myprintf("Skipping first %d lines\n",n_linestoskip,filename);
  for (i=0;i<n_linestoskip;i++)
    readsuccess=fgets(tline1,MAX_LINELENGTH,f1);


  /*------------------------------
    Read in the rest of the file
    ------------------------------*/
  chr_index=0;
  readstarts=(int *)(mappings[chr_index][index_readstart]);
  readends=(int *)(mappings[chr_index][index_readend]);
  readstrands=(signed char *)(mappings[chr_index][index_readstrand]);
  readscores=(float *)(mappings[chr_index][index_readscore]);
  linenumbers=(int *)(mappings[chr_index][index_linenumber]);
  nreads=*((int *)(mappings[chr_index][index_nreads]));
  maxreads=*((int *)(mappings[chr_index][index_maxreads]));

  nlines=0;
  readsuccess=fgets(tline1,MAX_LINELENGTH,f1);
  while (readsuccess != NULL)
  {
    nlines=nlines+1;
    if ((nlines%100000)==0)
    {
      myprintf("Reading line %d\n",nlines);
      myprintf("%s", tline1);
    }

    /*------------------------------
      parse the line, 6 columns of interest
      ------------------------------*/

    /* chromosome index 
    myprintf("Reading chromosome index\n"); */
    
    i1 = 0; 
    i2 = i1;
    while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 1 */
    memcpy(tempbuffer,tline1+i1,i2-i1);
    tempbuffer[i2-i1]='\0';
    /*myprintf("[%s]\n",tempbuffer);*/
    line_chrindex=-1;
    for (k=0;k<n_chromosomes;k++) 
      if (!strcmp(chromosomenames[k],tempbuffer))
      {
	line_chrindex=k;
	break;
      }
    if (line_chrindex < 0)
    {
      myprintf("Did not find a match for chromosome identifier [%s]\n", tempbuffer);
      return(mappings);
    }

    /*myprintf("Read chromosome index %d\n", line_chrindex);*/

    /* mapping start 
    myprintf("Reading mapping start\n");*/
    
    i1 = i2;
    while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 2 */
    i2 = i1;
    while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 2 */
    memcpy(tempbuffer,tline1+i1,i2-i1);
    tempbuffer[i2-i1]='\0';
    line_readstart = atoi(tempbuffer);
    /*myprintf("Read mapping start %d\n", line_readstart);*/

    /* mapping end 
    myprintf("Reading mapping end\n");*/
    
    i1 = i2;
    while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 3 */
    i2 = i1;
    while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 3 */
    memcpy(tempbuffer,tline1+i1,i2-i1);
    tempbuffer[i2-i1]='\0';
    line_readend = atoi(tempbuffer);
    /*myprintf("Read mapping end %d\n", line_readend);*/

    /* read name, we do not store this 
    myprintf("Reading mapping name\n");*/
    
    i1 = i2;
    while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 4 */
    i2 = i1;
    while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 4 */

    /* mapping score 
    myprintf("Reading mapping score\n");*/
    
    i1 = i2;
    while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 5 */
    i2 = i1;
    while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 5 */
    memcpy(tempbuffer,tline1+i1,i2-i1);
    tempbuffer[i2-i1]='\0';
    line_readscore = strtod(tempbuffer, NULL);
    /*myprintf("Read mapping score %f\n", line_readscore);*/

    /* strand identifier, '+' or '-' 
    myprintf("Reading mapping strand\n");*/
    
    i1 = i2;
    while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 6 */
    i2 = i1;
    while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 6 */
    if (tline1[i1] == '+') line_readstrand = 1; 
    else {if (tline1[i1] == '-') line_readstrand = -1; else line_readstrand = 0; }
    /*myprintf("Read mapping strand %d\n", line_readstrand);*/

    /* start working on the chromosome identified on the line, if we weren't already working on it */
    /*myprintf("Starting to work on correct chromosome\n");*/
    
    if (line_chrindex != chr_index)
    {
      /*myprintf("before: chr %d, nreads %d, maxreads %d",chr_index,nreads,maxreads);*/
      mappings[chr_index][index_readstart] = readstarts;
      mappings[chr_index][index_readend] = readends;
      mappings[chr_index][index_readstrand] = readstrands;
      mappings[chr_index][index_readscore] = readscores;
      mappings[chr_index][index_linenumber] = linenumbers;
      *((int *)(mappings[chr_index][index_nreads])) = nreads;
      *((int *)(mappings[chr_index][index_maxreads])) = maxreads;
    
      chr_index = line_chrindex;
      readstarts=(int *)(mappings[chr_index][index_readstart]);
      readends=(int *)(mappings[chr_index][index_readend]);
      readstrands=(signed char *)(mappings[chr_index][index_readstrand]);
      readscores=(float *)(mappings[chr_index][index_readscore]);
      linenumbers=(int *)(mappings[chr_index][index_linenumber]);
      nreads=*((int *)(mappings[chr_index][index_nreads]));
      maxreads=*((int *)(mappings[chr_index][index_maxreads]));
      /*myprintf("after: chr %d, nreads %d, maxreads %d",chr_index,nreads,maxreads);*/
    }

    /* add the mapped read to the data of the chromosome */
    /*myprintf("Adding the line to the data for the chromosome %d, %d reads, %d maxreads\n", chr_index, nreads,maxreads);*/
    
    nreads=nreads+1;    
    if (nreads>maxreads)
    {
      maxreads=nreads*2;
      readstarts=(int *)myrealloc(readstarts, sizeof(int)*maxreads);
      readends=(int *)myrealloc(readends, sizeof(int)*maxreads);
      readstrands=(signed char *)myrealloc(readstrands, sizeof(signed char)*maxreads);
      readscores=(float *)myrealloc(readscores, sizeof(float)*maxreads);
      linenumbers=(int *)myrealloc(linenumbers, sizeof(int)*maxreads);
    }
  
    readstarts[nreads-1]=line_readstart;
    readends[nreads-1]=line_readend;
    readstrands[nreads-1]=line_readstrand;
    readscores[nreads-1]=line_readscore;
    linenumbers[nreads-1]=nlines;  /* line number where this read occurred in the file */

    /* read in the next line */
    /*myprintf("Reading next line\n");*/
    readsuccess=fgets(tline1,MAX_LINELENGTH,f1);

    /* debug early stopping point */
    /* if (nlines > 1000) readsuccess = NULL; */
  }

  myprintf("Finished reading file, storing working data\n", tempbuffer);
  

  /* store the data of the last chromosome we were working on */
  mappings[chr_index][index_readstart] = readstarts;
  mappings[chr_index][index_readend] = readends;
  mappings[chr_index][index_readstrand] = readstrands;
  mappings[chr_index][index_readscore] = readscores;
  mappings[chr_index][index_linenumber] = linenumbers;
  *((int *)(mappings[chr_index][index_nreads])) = nreads;
  *((int *)(mappings[chr_index][index_maxreads])) = maxreads;

  myprintf("Discarding empty buffer data\n", tempbuffer);


  /* for all chromosomes, discard empty data */
  for (chr_index=0; chr_index<n_chromosomes;chr_index++)
  {
    readstarts=mappings[chr_index][index_readstart];
    readends=mappings[chr_index][index_readend];
    readstrands=mappings[chr_index][index_readstrand];
    readscores=mappings[chr_index][index_readscore];
    linenumbers=mappings[chr_index][index_linenumber];
    nreads=*((int *)(mappings[chr_index][index_nreads]));
    maxreads=*((int *)(mappings[chr_index][index_maxreads]));

    maxreads=nreads;
    readstarts=(int *)myrealloc(readstarts, sizeof(int)*maxreads);
    readends=(int *)myrealloc(readends, sizeof(int)*maxreads);
    readstrands=(signed char *)realloc(readstrands, sizeof(signed char)*maxreads);
    readscores=(float *)myrealloc(readscores, sizeof(float)*maxreads);
    linenumbers=(int *)myrealloc(linenumbers, sizeof(int)*maxreads);

    mappings[chr_index][index_readstart] = readstarts;
    mappings[chr_index][index_readend] = readends;
    mappings[chr_index][index_readstrand] = readstrands;
    mappings[chr_index][index_readscore] = readscores;
    mappings[chr_index][index_linenumber] = linenumbers;
    *((int *)(mappings[chr_index][index_nreads])) = nreads;
    *((int *)(mappings[chr_index][index_maxreads])) = maxreads;
  }

  myprintf("Closing file\n", tempbuffer);

  fclose(f1);

  return(mappings);
}


/*
int main(int argc, char *argv[])
{
  void ***mappings;
  mappings = read_mappingfile(argv[1]);
}
*/




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *tempfilename;
  char filename[MAX_LINELENGTH];

  mxArray *tempmappings;
  mxArray *tempreadstarts;
  mxArray *tempreadends;
  mxArray *tempreadstrands;
  mxArray *tempreadscores;
  mxArray *templinenumbers;
  mxArray *tempnreads;
  mxArray *tempmaxreads;
  int *tempreadstarts_data;
  int *tempreadends_data;
  signed char *tempreadstrands_data;
  float *tempreadscores_data;
  int *templinenumbers_data;
  int *tempnreads_data;
  int *tempmaxreads_data;


  void ***mappings;
  int nreads;
  int maxreads;
  int *readstarts;
  int *readends;
  signed char *readstrands;
  float *readscores;
  int *linenumbers;

  int chr_index, j, k;
  int tempindex[2];

  tempfilename = prhs[0];
  mxGetString(tempfilename, filename, MAX_LINELENGTH);

  myprintf("Trying to read from file [%s]\n",filename);
  mappings = read_mappingfile(filename);

  tempmappings = mxCreateCellMatrix(n_chromosomes,7);

  for (chr_index=0;chr_index<n_chromosomes;chr_index++)
  {    
    readstarts=mappings[chr_index][index_readstart];
    readends=mappings[chr_index][index_readend];
    readstrands=mappings[chr_index][index_readstrand];
    readscores=mappings[chr_index][index_readscore];
    linenumbers=mappings[chr_index][index_linenumber];
    nreads=*((int *)(mappings[chr_index][index_nreads]));
    maxreads=*((int *)(mappings[chr_index][index_maxreads]));

    tempreadstarts = mxCreateNumericMatrix(nreads,1,mxINT32_CLASS,0);
    tempreadstarts_data = (int *)mxGetData(tempreadstarts);

    tempreadends = mxCreateNumericMatrix(nreads,1,mxINT32_CLASS,0);
    tempreadends_data = (int *)mxGetData(tempreadends);

    tempreadstrands = mxCreateNumericMatrix(nreads,1,mxINT8_CLASS,0);
    tempreadstrands_data = (signed char *)mxGetData(tempreadstrands);

    tempreadscores = mxCreateNumericMatrix(nreads,1,mxSINGLE_CLASS,0);
    tempreadscores_data = (float *)mxGetData(tempreadscores);

    templinenumbers = mxCreateNumericMatrix(nreads,1,mxINT32_CLASS,0);
    templinenumbers_data = (int *)mxGetData(templinenumbers);

    tempnreads = mxCreateNumericMatrix(1,1,mxINT32_CLASS,0);
    tempnreads_data = (int *)mxGetData(tempnreads);

    tempmaxreads = mxCreateNumericMatrix(1,1,mxINT32_CLASS,0);
    tempmaxreads_data = (int *)mxGetData(tempmaxreads);

    tempnreads_data[0] = nreads;
    tempmaxreads_data[0] = maxreads;
    for (j=nreads-1;j>=0;j--)
    {
      tempreadstarts_data[j] = readstarts[j];
      tempreadends_data[j] = readends[j];
      tempreadstrands_data[j] = readstrands[j];
      tempreadscores_data[j] = readscores[j];
      templinenumbers_data[j] = linenumbers[j];
    }

    tempindex[0]=chr_index;
    tempindex[1]=index_readstart;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempreadstarts);

    tempindex[0]=chr_index;
    tempindex[1]=index_readend;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempreadends);

    tempindex[0]=chr_index;
    tempindex[1]=index_readstrand;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempreadstrands);

    tempindex[0]=chr_index;
    tempindex[1]=index_readscore;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempreadscores);

    tempindex[0]=chr_index;
    tempindex[1]=index_nreads;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempnreads);

    tempindex[0]=chr_index;
    tempindex[1]=index_maxreads;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempmaxreads);

    tempindex[0]=chr_index;
    tempindex[1]=index_linenumber;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,templinenumbers);

    mxFree(readstarts);
    mxFree(readends);
    mxFree(readscores);
    mxFree(linenumbers);
  }

  plhs[0] = tempmappings;
  
}


