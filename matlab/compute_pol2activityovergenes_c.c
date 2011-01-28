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



void compute_overlaps(int n_reads,int *tempreadstarts_data,int *tempreadends_data,signed char *tempreadstrands_data,float *tempreadscores_data,int fragment_length,int n_allowed_duplicates,int n_bins,int *tempbinstarts_data,int *tempbinends_data,signed char *tempbinstrands_data,int *n_subbins, double **subbins_data, int subbin_length)
{
  int readstart, readend, tempreadstart, tempreadend, modified_readstart, modified_readend;
  int tempstart,tempend;
  signed char readstrand;
  float readscore;
  double binscoredelta;
  int n_duplicatesfound;
  int i, j, k, k0;
  int l, m, firstsubbin, lastsubbin;
  double firstsubbin_proportion, lastsubbin_proportion;
  int verbose=0;

  readstart=-1;
  readend=-1;
  n_duplicatesfound=0;

  for (i=0;i<n_reads;i++)
  {
    verbose=0;
    if ((i%50000)==0) verbose=1;

    tempreadstart=tempreadstarts_data[i];
    tempreadend=tempreadends_data[i];


    if (verbose==1)
      myprintf("Read %d of %d, start %d, end %d\n", i, n_reads, tempreadstart, tempreadend);

  
    /*-------------------------------------
      check for a duplicate (assume reads are sorted so duplicates
      are listed one after another)
      -------------------------------------*/
    if ((tempreadstart==readstart) && (tempreadend==readend))
    {
      n_duplicatesfound=n_duplicatesfound+1;
      if (verbose==1)
	myprintf("Found duplicate no. %d of this read\n", n_duplicatesfound);
    }
    else
    {
      n_duplicatesfound=1;
      readstart=tempreadstart;
      readend=tempreadend;
      if (verbose==1)
	myprintf("First instance of this read\n", n_duplicatesfound);
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

      if (verbose==1)
	myprintf("Extended location: %d - %d\n", modified_readstart, modified_readend);
    

      /*-------------------------------------
	% Compute assignment to bins.
	% Assume the bins have been sorted according to start index.
	% Assume that bins may overlap. A read that falls into multiple bins
	% is currently counted separately for each bin (i.e. the
	% influence of the read is multiplied).
	%-------------------------------------*/

      /* find latest bin (greatest index) whose start is not after the read */
      k0=n_bins-1; while ((k0>=0) && (tempbinstarts_data[k0]>modified_readend)) k0=k0-1;

      if (verbose==1)
	myprintf("Extended location: %d - %d\n", modified_readstart, modified_readend);
    
      for (k=k0;k>=0;k--)
      {
	/*
	  % Since we know the end of the read is at or after the start
	  % of the bin, all we need to check is that the start of the 
	  % read has not happened after the bin (i.e. check that the
	  % whole read is not to the right of the bin).
	*/

	/*
        if (verbose==1)
	    myprintf("Read %d of %d (%d - %d): testing bin %d of %d\n", 
		     i, n_reads, tempreadstart, tempreadend, k, n_bins);
	*/


	if (modified_readstart <= tempbinends_data[k])
	{
	  if (verbose==1)
	    myprintf("Read %d of %d (%d - %d) matched bin %d of %d (%d - %d)\n", 
		     i, n_reads, tempreadstart, tempreadend, k, n_bins, tempbinstarts_data[k], tempbinends_data[k]);

      
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

	  if (1)
	  {	        
	    /*-------------------------------------
	      % Assign the read to subbins inside the bin.
	      % Currently, amount added to each subbin depends on length of
	      % overlap and on score of the read.
              % Note that orientation of the subbins depends on the strand of
              % the bin! We especially want to model the *end* of the bin, so 
              % subbin orientations are reversed: 
              %     If strand=1, the first subbin is at the end of the bin.
              %     If strand=-1, the first subbin is at the start of the bin.
	      -------------------------------------*/
	    binscoredelta=readscore;

	    if (tempbinstrands_data[k] < 0)
	    {
	      /* Count subbins from the start of the bin */
	      firstsubbin=(tempstart-tempbinstarts_data[k])/subbin_length;
	      lastsubbin=(tempend-tempbinstarts_data[k])/subbin_length;
	    }
	    else if (tempbinstrands_data[k] > 0)
	    {
	      /* Count subbins from the end of the bin */
	      firstsubbin=(tempbinends_data[k]-tempend)/subbin_length;
	      lastsubbin=(tempbinends_data[k]-tempstart)/subbin_length;
	    }

	    if ((firstsubbin<0)||(lastsubbin<0)||(firstsubbin>=n_subbins[k])||(lastsubbin>=n_subbins[k]))
	    {
	      myprintf("PROBLEM: Read %d of %d (%d - %d) matched bin %d of %d (%d - %d), subbins %d-%d of %d\n", 
		       i, n_reads, tempreadstart, tempreadend, k, n_bins, tempbinstarts_data[k], tempbinends_data[k], firstsubbin, lastsubbin, n_subbins[k]);
	      return;
	    }

	    if (verbose==1)
	      myprintf("Read %d of %d (%d - %d) matched bin %d of %d (%d - %d), subbins %d-%d\n", 
		       i, n_reads, tempreadstart, tempreadend, k, n_bins, tempbinstarts_data[k], tempbinends_data[k], firstsubbin, lastsubbin);
	  }

	  if (1)
	  {
	    if (lastsubbin > firstsubbin)
	    {
	      if (tempbinstrands_data[k] < 0)
	      {
		/* Count subbins from the start of the bin */
		firstsubbin_proportion=tempbinstarts_data[k]+(firstsubbin+1)*subbin_length - tempstart;
		lastsubbin_proportion=tempend - tempbinstarts_data[k] - lastsubbin*subbin_length+1;
	      }
	      else if (tempbinstrands_data[k] > 0)
	      {
		/* Count subbins from the end of the bin */
		firstsubbin_proportion=tempend - tempbinends_data[k] + (firstsubbin+1)*subbin_length;
		lastsubbin_proportion=tempbinends_data[k] - lastsubbin*subbin_length - tempstart+1;
	      }

	      subbins_data[k][firstsubbin] += firstsubbin_proportion*binscoredelta;
	      subbins_data[k][lastsubbin] += lastsubbin_proportion*binscoredelta;
	      for (l=firstsubbin+1;l<lastsubbin;l++)
		subbins_data[k][l] += subbin_length*binscoredelta;	    

	      if ((firstsubbin<0) || (lastsubbin>n_subbins[k]) || (firstsubbin_proportion<0) || (lastsubbin_proportion<0))
		myprintf("PROBLEM: Read %d of %d (%d - %d) matched bin %d of %d (%d - %d), subbins %d-%d, firstprop %f, lastprop %f\n", 
			 i, n_reads, tempreadstart, tempreadend, k, n_bins, tempbinstarts_data[k], tempbinends_data[k], firstsubbin, lastsubbin, firstsubbin_proportion, lastsubbin_proportion);
	    }
	    else
	    {
	      subbins_data[k][firstsubbin] += (tempend-tempstart+1)*binscoredelta;
	    }
	  }
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
  mxArray *tempsubbinlength;

  mxArray *tempnbins;
  mxArray *tempbinstarts;
  mxArray *tempbinends;
  mxArray *tempbinstrands;
  mxArray *tempsubbins;
  mxArray **tempbinheights;

  mxArray *tempchrindex;
  mxArray *tempfragmentlength;
  mxArray *tempnallowedduplicates;

  int *tempreadstarts_data;
  int *tempreadends_data;
  signed char *tempreadstrands_data;
  float *tempreadscores_data;
  int *tempnreads_data;
  int *tempmaxreads_data;
  int *tempsubbinlength_data;

  int *tempnbins_data;
  int *tempbinstarts_data;
  int *tempbinends_data;
  signed char *tempbinstrands_data;
  int *n_subbins;
  double **subbins_data;

  int *tempchrindex_data;
  int *tempfragmentlength_data;
  int *tempnallowedduplicates_data;
  int tempindex[2];

  int chr_index, n_reads, fragment_length, n_allowed_duplicates, n_bins, subbin_length;
  int j, k;

  tempmappings = prhs[0];
  tempchrindex = prhs[1];
  tempnbins = prhs[2];
  tempbinstarts = prhs[3];
  tempbinends = prhs[4];
  tempbinstrands = prhs[5];
  tempfragmentlength = prhs[6];
  tempnallowedduplicates = prhs[7];
  tempsubbinlength = prhs[8];

  /* Get chromosome index, convert to zero-based index */
  tempchrindex_data = (int *)mxGetData(tempchrindex);
  chr_index = tempchrindex_data[0];
  chr_index=chr_index-1;

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
  tempreadscores_data = (float *)mxGetData(tempreadscores);

  /* Get number of bins and vectors of bin starts and bin ends and bin strands */
  tempnbins_data = (int *)mxGetData(tempnbins);
  n_bins=tempnbins_data[0];
  tempbinstarts_data = (int *)mxGetData(tempbinstarts);
  tempbinends_data = (int *)mxGetData(tempbinends);
  tempbinstrands_data = (signed char *)mxGetData(tempbinstrands);

  if (1)
  {
    /* Get desired length of sub-bins */
    tempsubbinlength_data = (int *)mxGetData(tempsubbinlength);
    subbin_length=tempsubbinlength_data[0];

    /* Create output cell array for subbin heights */
    tempsubbins = mxCreateCellMatrix(n_bins,1);
    n_subbins = (int *)mymalloc(n_bins*sizeof(int));
    subbins_data = (double **)mymalloc(n_bins*sizeof(double *));
    tempbinheights = (mxArray **)mymalloc(n_bins*sizeof(mxArray *));
    for (j=0; j<n_bins; j++)
    {
      /* Compute number of subbins, rounding up */
      k = (tempbinends_data[j]-tempbinstarts_data[j]+1)/subbin_length; 
      if (tempbinends_data[j]-tempbinstarts_data[j]+1 > k*subbin_length)
	k=k+1;
      tempbinheights[j] = mxCreateNumericMatrix(k,1,mxDOUBLE_CLASS,0);
      subbins_data[j] = (double *)mxGetData(tempbinheights[j]);
      n_subbins[j] = k;
    }
  }
  /* Compute overlaps */

  
  myprintf("%d reads %p %p %p %p, fragmentlength %d, duplicates %d, bins %d %p %p, subbins %p %p length %d\n",
	   n_reads,tempreadstarts_data,tempreadends_data,tempreadstrands_data,tempreadscores_data,fragment_length,n_allowed_duplicates,n_bins,tempbinstarts_data,tempbinends_data,n_subbins,subbins_data,subbin_length);

  compute_overlaps(n_reads,tempreadstarts_data,tempreadends_data,tempreadstrands_data,tempreadscores_data,fragment_length,n_allowed_duplicates,n_bins,tempbinstarts_data,tempbinends_data,tempbinstrands_data,n_subbins,subbins_data,subbin_length);
  

  if (1)
  {
    /* Place answers into cell array */
    for (j=0; j<n_bins; j++)
    {
      tempindex[0]=j;
      tempindex[1]=0;
      k = mxCalcSingleSubscript(tempsubbins,2,tempindex);
      mxSetCell(tempsubbins,k,tempbinheights[j]);
    }
    plhs[0] = tempsubbins;

    /* Clean up */
    mxFree(n_subbins);
    mxFree(subbins_data);  
    mxFree(tempbinheights);
  }

}


