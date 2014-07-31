#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

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




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *polbins;
  mxArray *rnabins;
  mxArray *nbins;
  mxArray *ntimepoints;
  mxArray *corrcutoff;
  mxArray *ntrials;
  double *polbins_data;
  double *rnabins_data;
  int *nbins_data;
  int *ntimepoints_data;
  int *ntrials_data;
  double *corrcutoff_data;

  mxArray *trialresults;
  int *trialresults_data;


  int n_bins, n_timepoints, n_trials;
  double corr_cutoff;

  int trial;
  int i, j, k, l, m, n;
  int tempindex[2];

  int *myperm1, *myperm2;
  double *polvariances, *rnavariances;
  double *polmeans, *rnameans;
  int tempresult;
  double tempmean, tempvar, tempvalue, tempcorr;


  polbins = prhs[0];
  rnabins = prhs[1];
  nbins = prhs[2];
  ntimepoints = prhs[3];
  corrcutoff = prhs[4];
  ntrials = prhs[5];


  /* Get POL and RNA bins and other inputs from MATLAB */
  polbins_data = (double *)mxGetData(polbins);
  rnabins_data = (double *)mxGetData(rnabins);
  nbins_data = (int *)mxGetData(nbins);
  n_bins = nbins_data[0];
  ntimepoints_data = (int *)mxGetData(ntimepoints);
  n_timepoints = ntimepoints_data[0];
  corrcutoff_data = (double *)mxGetData(corrcutoff);
  corr_cutoff = corrcutoff_data[0];
  ntrials_data = (int *)mxGetData(ntrials);
  n_trials = ntrials_data[0];

  /* Create output vector of trial results */
  trialresults = mxCreateNumericMatrix(n_bins,1,mxINT32_CLASS,0);
  trial_results = (double *)mxGetData(tempbinheights);
  plhs[0] = trialresults;



  /* Matlab uses column-major order!!! A(i,j) --> A[j*n_i + i]  */
#define BINELEMENT(myrow,mycol) (mycol*n_bins+myrow)

  /* compute means and (squared roots of) variances for pol and rna */
  polvariances=(double *)mymalloc(n_bins*sizeof(double));
  rnavariances=(double *)mymalloc(n_bins*sizeof(double));
  polmeans=(double *)mymalloc(n_bins*sizeof(double));
  rnameans=(double *)mymalloc(n_bins*sizeof(double));
  for (i=0;i<n_bins;i++)
  {
    tempmean=0;
    for (j=0;j<n_timepoints;j++)
      tempmean+=polbins_data[BINELEMENT(i,j)];
    polmeans[i]=tempmean/n_timepoints;

    tempmean=0;
    for (j=0;j<n_timepoints;j++)
      tempmean+=rnabins_data[BINELEMENT(i,j)];
    rnameans[i]=tempmean/n_timepoints;

    tempmean=polmeans[i];
    tempvar=0;
    for (j=0;j<n_timepoints;j++)
    {
      tempvalue=polbins_data[BINELEMENT(i,j)]-tempmean;
      tempvar+=tempvalue*tempvalue;
    }
    polvars[i]=tempvar/n_timepoints;
    polvars[i]=sqrt(polvars[i]);

    tempmean=rnameans[i];
    tempvar=0;
    for (j=0;j<n_timepoints;j++)
    {
      tempvalue=rnabins_data[BINELEMENT(i,j)]-tempmean;
      tempvar+=tempvalue*tempvalue;
    }
    rnavars[i]=tempvar/n_timepoints;
    rnavars[i]=sqrt(rnavars[i]);
  }

  trialresults_data = (int *)mymalloc(n_trials*sizeof(int)); 
  myperm1=(int *)mymalloc(n_bins*sizeof(int)); 
  myperm2=(int *)mymalloc(n_bins*sizeof(int)); 

  for (trial=0;trial<ntrials;trial++)
  {
    tempresult=0;
    if ((trial%100)==0) 
      myprintf("trial %d\n",trial);    

    /* Create random permutation */
    for (k = 0; k < n_bins; k++)
      myperm1[k]=k;
    for (k = 0; k < n_bins; k++)
    {
      l=random()%(n_bins-k);
      myperm2[k]=myperm1[l];
      m=myperm1[n_bins-k];
      myperm1[n_bins-k]=myperm1[l];
      myperm1[l]=m;
    }

    /* Compute correlations, check if they are above the cutoff */
    for (i=0;i<n_bins;i++)
    {
      tempcorr=0;
      k=myperm2[i];
      for (j=0;j<n_timepoints;j++)
	tempcorr+=(polbins_data[BINELEMENT(i,j)]-polmeans[i])*(rnabins_data[BINELEMENT(k,j)]-rnameans[k]);
      tempcorr/=polvars[i]*rnavars[k];
      if (tempcorr>=corr_cutoff) tempresult++;
    }

    trialresults_data[trial]=tempresult;
  }

  myfree(myperm1);
  myfree(myperm2);
  myfree(rnameans);
  myfree(rnavars);
  myfree(polmeans);
  myfree(polvars);  
}


