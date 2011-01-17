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

#define n_chromosomes 25
const char *chromosomenames[]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT"};

#define index_readstart 0
#define index_readend 1
#define index_readstrand 2
#define index_readscore 3
#define index_nreads 4
#define index_maxreads 5
#define index_linenumber 6




/*
int main(int argc, char *argv[])
{
  void ***mappings;
  mappings = read_mappingfile(argv[1]);
}
*/



void compute_overlaps(int n_reads,int *tempreadstarts_data,int *tempreadends_data,signed char *tempreadstrands_data,float *tempreadscores_data,int fragment_length,int n_allowed_duplicates,int n_bins,int *tempbinstarts_data,int *tempbinends_data,double *tempbinheights_data)
{
  int readstart, readend, tempreadstart, tempreadend, modified_readstart, modified_readend;
  int tempstart,tempend;
  signed char readstrand;
  float readscore;
  int n_duplicatesfound;
  int i, j, k, k0;


  readstart=-1;
  readend=-1;
  n_duplicatesfound=0;

  for (i=0;i<n_reads;i++)
  {
    if ((i%50000)==0)
      myprintf("Read %d\n", i);
  
    tempreadstart=tempreadstarts_data[i];
    tempreadend=tempreadends_data[i];
  
    /*-------------------------------------
      check for a duplicate (assume reads are sorted so duplicates
      are listed one after another)
      -------------------------------------*/
    if ((tempreadstart==readstart) && (tempreadend==readend))
    {
      n_duplicatesfound=n_duplicatesfound+1;
    }
    else
    {
      n_duplicatesfound=1;
      readstart=tempreadstart;
      readend=tempreadend;
    }
  
    if (n_duplicatesfound <= n_allowed_duplicates)
    {
      readstrand=tempreadstrands_data[i];
      readscore=tempreadscores_data[i];

      /*-------------------------------------
	shift and extend the read using the desired fragment length
	-------------------------------------*/
      if (fragment_length > 0)
      {  
	if (readstrand>0)
	{
	  /* OLD VERSION
	   shift 
	  modified_readstart=readstart+fragment_length/2;
	  modified_readend=readend+fragment_length/2;
  
	   extend 
	  modified_readstart=modified_readstart-fragment_length/2;
	  modified_readend=modified_readstart+fragment_length;  
	  */

	  /* NEW version, combines shift and extend */
	  modified_readstart=readstart;
	  modified_readend=readstart+fragment_length;
	}
      
	if (readstrand<0)
	{
	  /* OLD VERSION
	     % shift
	     modified_readstart=readstart-fragment_length/2;
	     modified_readend=readend-fragment_length/2;
	
	     % extend
	     modified_readend=modified_readend+fragment_length/2;  
	     modified_readstart=modified_readend-fragment_length;
	  */
	  
	  /* NEW version, combines shift and extend */
	  modified_readend=readend;
	  modified_readstart=readend-fragment_length;
	}
      }
    

      /*-------------------------------------
	% Compute assignment to bins.
	% Assume the bins have been sorted according to start index.
	% Assume that bins may overlap. A read that falls into multiple bins
	% is currently counted separately for each bin (i.e. the
	% influence of the read is multiplied).
	%-------------------------------------*/

      /* find latest bin (greatest index) whose start is not after the read */
      k0=n_bins-1; while ((k0>=0) && (tempbinstarts_data[k0]>modified_readend)) k0=k0-1;
    
      for (k=k0;k>=0;k--)
      {
	/*
	  % Since we know the end of the read is at or after the start
	  % of the bin, all we need to check is that the start of the 
	  % read has not happened after the bin (i.e. check that the
	  % whole read is not to the right of the bin).
	*/
	if (modified_readstart <= tempbinends_data[k])
	{
      
	  /*-------------------------------------
	    % Compute overlapping portion between the read and the bin.
	    -------------------------------------*/
	  if (modified_readend>=tempbinends_data[k])
	    tempend=tempbinends_data[k]; 
	  else 
	    tempend=modified_readend;
	 
	  if (modified_readstart<=tempbinstarts_data[k])
	    tempstart=tempbinstarts_data[k];
	  else
	    tempstart=modified_readstart;      

	  if (tempstart>tempend)
	    myprintf("problem at read %d, bin %d: bin %d, %d, read %d, %d, overlap %d\n",
		     i, k, tempbinstarts_data[k],tempbinends_data[k],modified_readstart,
		modified_readend,tempend-tempstart+1);
	        
	  /*-------------------------------------
	    % Assign the read to the bin.
	    % Currently, amount added to the bin depends on length of
	    % overlap and on score of the read.      
	    -------------------------------------*/
	  tempbinheights_data[k]=tempbinheights_data[k]+(tempend-tempstart+1)*readscore;
	}
      }
    }
    
  }

  
}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *tempfilename1;
  mxArray *tempfilename2;
  char filename1[MAX_LINELENGTH];
  char filename2[MAX_LINELENGTH];

  mxArray *tempmappings;
  mxArray *tempreadstarts;
  mxArray *tempreadends;
  mxArray *tempreadstrands;
  mxArray *tempreadscores;
  mxArray *tempnreads;
  mxArray *tempmaxreads;

  mxArray *tempnbins;
  mxArray *tempbinstarts;
  mxArray *tempbinends;
  mxArray *tempbinheights;

  mxArray *tempchrindex;
  mxArray *tempfragmentlength;
  mxArray *tempnallowedduplicates;

  int *tempreadstarts_data;
  int *tempreadends_data;
  signed char *tempreadstrands_data;
  float *tempreadscores_data;
  int *tempnreads_data;
  int *tempmaxreads_data;

  int *tempnbins_data;
  int *tempbinstarts_data;
  int *tempbinends_data;
  double *tempbinheights_data;

  int *tempchrindex_data;
  int *tempfragmentlength_data;
  int *tempnallowedduplicates_data;

  int chr_index, n_reads, fragment_length, n_allowed_duplicates, n_bins;
  int j, k;
  int tempindex[2];



  tempmappings = prhs[0];
  tempchrindex = prhs[1];
  tempnbins = prhs[2];
  tempbinstarts = prhs[3];
  tempbinends = prhs[4];
  tempfragmentlength = prhs[5];
  tempnallowedduplicates = prhs[6];

  /* Get chromosome index */
  tempchrindex_data = (int *)mxGetData(tempchrindex);
  chr_index = tempchrindex_data[0];

  /* Get desired fragment length */
  tempfragmentlength_data = (int *)mxGetData(tempfragmentlength);
  fragment_length = tempfragmentlength_data[0];

  /* Get maximum number of allowed duplicates */
  tempnallowedduplicates_data = (int *)mxGetData(tempnallowedduplicates);
  n_allowed_duplicates = tempnallowedduplicates_data[0];

  /* Get number of reads */
  tempindex[0]=chr_index;
  tempindex[1]=index_nreads;
  j = mxCalcSingleSubscript(tempmappings,2,tempindex);
  tempnreads=mxGetCell(tempmappings,j);
  tempnreads_data = (int *)mxGetData(tempnreads);
  n_reads = tempnreads_data[0];

  /* Get vector of read starts */
  tempindex[0]=chr_index;
  tempindex[1]=index_readstart;
  j = mxCalcSingleSubscript(tempmappings,2,tempindex);
  tempreadstarts=mxGetCell(tempmappings,j);
  tempreadstarts_data = (int *)mxGetData(tempreadstarts);

  /* Get vector of read ends */
  tempindex[0]=chr_index;
  tempindex[1]=index_readend;
  j = mxCalcSingleSubscript(tempmappings,2,tempindex);
  tempreadends=mxGetCell(tempmappings,j);
  tempreadends_data = (int *)mxGetData(tempreadends);

  /* Get vector of read strands */
  tempindex[0]=chr_index;
  tempindex[1]=index_readstrand;
  j = mxCalcSingleSubscript(tempmappings,2,tempindex);
  tempreadstrands=mxGetCell(tempmappings,j);
  tempreadstrands_data = (signed char *)mxGetData(tempreadstrands);

  /* Get vector of read scores */
  tempindex[0]=chr_index;
  tempindex[1]=index_readscore;
  j = mxCalcSingleSubscript(tempmappings,2,tempindex);
  tempreadscores=mxGetCell(tempmappings,j);
  tempreadscores_data = (int *)mxGetData(tempreadscores);

  /* Get number of bins and vectors of bin starts and bin ends */
  tempnbins_data = (int *)mxGetData(tempnbins);
  n_bins=tempnbins_data[0];
  tempbinstarts_data = (int *)mxGetData(tempbinstarts);
  tempbinends_data = (int *)mxGetData(tempbinends);

  /* Create output vector of bin heights */
  tempbinheights = mxCreateNumericMatrix(n_bins,1,mxDOUBLE_CLASS,0);
  tempbinheights_data = (double *)mxGetData(tempbinheights);
  plhs[0] = tempbinheights;


  /* Compute overlaps */
  compute_overlaps(n_reads,tempreadstarts_data,tempreadends_data,tempreadstrands_data,tempreadscores_data,fragment_length,n_allowed_duplicates,n_bins,tempbinstarts_data,tempbinends_data,tempbinheights_data);


  
}


