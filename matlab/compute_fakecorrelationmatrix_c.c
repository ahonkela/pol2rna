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


/*
int main(int argc, char *argv[])
{
  void ***mappings;
  mappings = read_mappingfile(argv[1]);
}
*/





const double mytimepoints[10]={0, 5, 10, 20, 40, 80, 160, 320, 640, 1280};


void compute_fakecorrelationmatrix
(
  int n_genes, 
  int n_timepoints,
  int profilelength, 
  double ***allgenebins_data,
  double pol_spatialspeed,

  double *binmeans,
  double *binvariances,
  double *binsamples,
  double *corrmatrix,
  double *ncorrsamples
)
{
  int i, j, g, timeindex, timeindex2;
  double t, t2;
  double value_t, value_t2;

  myprintf("Processing %d genes, %d timepoints, profilelength %d, pol-spatialspeed %f\n", n_genes,n_timepoints,profilelength,pol_spatialspeed);


  /* compute bin means and variances */
  for (i=0;i<profilelength+1;i++)
  {
    binmeans[i]=0;
    binvariances[i]=0;
    binsamples[i]=0;
  }

  for (i=0;i<profilelength;i++) /* omit RNA for now */  
  {
    for (g=0;g<n_genes; g++)
    {
      for (timeindex=0;timeindex<n_timepoints;timeindex++)
      {
	value_t=allgenebins_data[g][timeindex][i];
	binmeans[i]=binmeans[i]+value_t;
	binsamples[i]=binsamples[i]+1;
      }
    }
  }
  for (i=0;i<profilelength;i++) /* omit RNA for now */  
  {
    binmeans[i]=binmeans[i]/binsamples[i];
  }

  /* compute bin variances */
  for (i=0;i<profilelength;i++) /* omit RNA for now */  
  {    
    for (g=0;g<n_genes; g++)
    {
      for (timeindex=0;timeindex<n_timepoints;timeindex++)
      {
	value_t=allgenebins_data[g][timeindex][i];
	binvariances[i]=binvariances[i]+(value_t-binmeans[i])*(value_t-binmeans[i]);
      }
    }
  }
  for (i=0;i<profilelength;i++) /* omit RNA for now */  
  {
    binvariances[i]=binvariances[i]/binsamples[i];
  }

  /* spatial speed of POL2 (with pauses and all) in bins/minute = (basepairs/minute)/binlength */
  pol_spatialspeed=(2000.0/5.0)/200.0;


#define ACCESSELEMENT(i,j) ((j)*(profilelength+1)+(i))
  for (i=0;i<profilelength+1;i++)
  {
    for (j=0;j<profilelength+1;j++)
    {
      corrmatrix[ACCESSELEMENT(i,j)]=0;
      ncorrsamples[ACCESSELEMENT(i,j)]=0;
    }
  }


  for (timeindex=0;timeindex<n_timepoints;timeindex++)
  {
    t=mytimepoints[timeindex];
    for (g=0;g<n_genes; g++)
    {
      myprintf("Processing timeindex %d (%f) of %d, gene %d of %d\n", timeindex, t, n_timepoints, g, n_genes);
      for (i=0;i<profilelength;i++) /* omit RNA for now */
      {
	value_t=allgenebins_data[g][timeindex][i];
            
	for (j=0;j<profilelength;j++) /* omit RNA for now */
        {
	  t2=t+(j-i)/pol_spatialspeed;
	  if ((t2 >= mytimepoints[0]) && (t2 <= mytimepoints[n_timepoints-1]))
	  {
	    timeindex2=0;
	    while((timeindex2<n_timepoints)&&(mytimepoints[timeindex2]<t2)) 
	      timeindex2++;
	    if (timeindex2 >= n_timepoints)
	    {
	      myprintf("ERROR: interpolated timepoint not found, t2=%f\n", t2);
	    }
	    
	    if (timeindex2<n_timepoints-1)
	    {
	      value_t2=allgenebins_data[g][timeindex2][j]+(t2-mytimepoints[timeindex2])/(mytimepoints[timeindex2+1]-mytimepoints[timeindex2])*(allgenebins_data[g][timeindex2+1][j]-allgenebins_data[g][timeindex2][j]);
	    }
	    else
	    {
	      value_t2=allgenebins_data[g][timeindex2][j];
	    }
	    corrmatrix[ACCESSELEMENT(i,j)]=corrmatrix[ACCESSELEMENT(i,j)]+value_t*value_t2;
	    ncorrsamples[ACCESSELEMENT(i,j)]=ncorrsamples[ACCESSELEMENT(i,j)]+1;
	  }
	} /* For all locations j */
      } /* For all locations i */
    } /* For all genes */
  } /* For all time indices */


  for (i=0;i<profilelength;i++) /* omit RNA for now */
  {
    for (j=0;j<profilelength;j++) /* omit RNA for now */
    {
      corrmatrix[ACCESSELEMENT(i,j)] /= ncorrsamples[ACCESSELEMENT(i,j)];
    }
  }
}





void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *tempfilename1;
  mxArray *tempfilename2;
  char filename1[MAX_LINELENGTH];
  char filename2[MAX_LINELENGTH];

  mxArray *tempngenes;
  int *tempngenes_data;
  int n_genes;

  mxArray *tempntimepoints;
  int *tempntimepoints_data;
  int n_timepoints;

  mxArray *tempprofilelength;
  int *tempprofilelength_data;
  int profilelength;

  mxArray *tempallgenebins;
  mxArray *tempprofile;
  double ***allgenebins_data;

  mxArray *temppol2spatialspeed;
  double *temppol2spatialspeed_data;
  double pol_spatialspeed;

  mxArray *tempbinmeans;
  double *binmeans_data;

  mxArray *tempbinvariances;
  double *binvariances_data;

  mxArray *tempbinsamples;
  double *binsamples_data;

  mxArray *tempcorrmatrix;
  double *corrmatrix_data;

  mxArray *tempcorrsamples;
  double *corrsamples_data;

  int tempindex[2];

  int i, j, k;



  tempngenes = prhs[0];
  tempngenes_data = (int *)mxGetData(tempngenes);
  n_genes = tempngenes_data[0];

  tempntimepoints = prhs[1];
  tempntimepoints_data = (int *)mxGetData(tempntimepoints);
  n_timepoints = tempntimepoints_data[0];

  tempprofilelength = prhs[2];
  tempprofilelength_data = (int *)mxGetData(tempprofilelength);
  profilelength = tempprofilelength_data[0];

  tempallgenebins = prhs[3];
  allgenebins_data=(double ***)mymalloc(n_genes*sizeof(double **));
  for (i=0;i<n_genes;i++)
  {
    allgenebins_data[i]=(double **)mymalloc(n_timepoints*sizeof(double *));
    for (j=0;j<n_timepoints;j++)
    {
      /* Get gene bin profile */
      tempindex[0]=i;
      tempindex[1]=j;
      k = mxCalcSingleSubscript(tempallgenebins,2,tempindex);
      tempprofile=mxGetCell(tempallgenebins,k);
      allgenebins_data[i][j]=(double *)mxGetData(tempprofile);
    }
  }

  temppol2spatialspeed = prhs[4];
  temppol2spatialspeed_data = (double *)mxGetData(temppol2spatialspeed);
  pol_spatialspeed = temppol2spatialspeed_data[0];

  
  tempbinmeans = mxCreateNumericMatrix(profilelength+1,1,mxDOUBLE_CLASS,0);
  binmeans_data = (double *)mxGetData(tempbinmeans);
  plhs[0]=tempbinmeans;

  tempbinvariances = mxCreateNumericMatrix(profilelength+1,1,mxDOUBLE_CLASS,0);
  binvariances_data = (double *)mxGetData(tempbinvariances);
  plhs[1]=tempbinvariances;

  tempbinsamples = mxCreateNumericMatrix(profilelength+1,1,mxDOUBLE_CLASS,0);
  binsamples_data = (double *)mxGetData(tempbinsamples);
  plhs[2]=tempbinsamples;

  tempcorrmatrix = mxCreateNumericMatrix(profilelength+1,profilelength+1,mxDOUBLE_CLASS,0);
  corrmatrix_data = (double *)mxGetData(tempcorrmatrix);
  plhs[3]=tempcorrmatrix;

  tempcorrsamples = mxCreateNumericMatrix(profilelength+1,profilelength+1,mxDOUBLE_CLASS,0);
  corrsamples_data = (double *)mxGetData(tempcorrsamples);
  plhs[4]=tempcorrsamples;

  compute_fakecorrelationmatrix(n_genes,n_timepoints,profilelength,allgenebins_data,pol_spatialspeed,binmeans_data,binvariances_data,binsamples_data,corrmatrix_data,corrsamples_data);

  /* Clean up */
  for (i=0;i<n_genes;i++)
  {
    myfree(allgenebins_data[i]);
  }
  myfree(allgenebins_data);
}


