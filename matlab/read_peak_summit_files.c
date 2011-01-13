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
#define n_linestoskip_peaks 0
#define n_linestoskip_summits 0

#define n_chromosomes 25
const char *chromosomenames[]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chrMT"};


#define index_peakstart 0
#define index_peakend 1
#define index_peakscore 2
#define index_summitstart 3
#define index_summitend 4
#define index_summitheight 5
#define index_linenumber 6
#define index_npeaks 7
#define index_maxpeaks 8


void ***read_peak_summit_files(char *peakfilename, char *summitfilename)
{
  FILE *f1, *f2;
  int peaks_available = 0;
  int summits_available = 0;
  char tline1[MAX_LINELENGTH];
  char tline2[MAX_LINELENGTH];
  char tempbuffer1[MAX_LINELENGTH];
  char tempbuffer2[MAX_LINELENGTH];
  char peakname1[MAX_LINELENGTH];
  char peakname2[MAX_LINELENGTH];
  char *readsuccess1, *readsuccess2;

  int i, k, i1, i2;
  int chr_index;
  int nlines;

  void ***mappings;
  int npeaks;
  int maxpeaks;
  int *peakstarts;
  int *peakends;
  float *peakscores;
  int *summitstarts;
  int *summitends;
  float *summitheights;
  int *linenumbers;

  int line_chrindex1;
  int line_chrindex2;
  int line_peakstart;
  int line_peakend;
  float line_peakscore;
  int line_summitstart;
  int line_summitend;
  float line_summitheight;


  /*------------------------------
  Initialize data cell array
  ------------------------------*/
  myprintf("Initializing cell array\n");

  mappings=(void ***)mymalloc(n_chromosomes*sizeof(void **));
  for (chr_index=0;chr_index<n_chromosomes;chr_index++)
  {
    mappings[chr_index] = (void **)mymalloc(9*sizeof(void *));

    mappings[chr_index][index_peakstart] = NULL;
    mappings[chr_index][index_peakend] = NULL;
    mappings[chr_index][index_peakscore] = NULL;
    mappings[chr_index][index_summitstart] = NULL;
    mappings[chr_index][index_summitend] = NULL;
    mappings[chr_index][index_summitheight] = NULL;
    mappings[chr_index][index_linenumber] = NULL;
    mappings[chr_index][index_npeaks] = mymalloc(1*sizeof(int));
    mappings[chr_index][index_maxpeaks] = mymalloc(1*sizeof(int));
    *((int *)(mappings[chr_index][index_npeaks])) = 0;
    *((int *)(mappings[chr_index][index_maxpeaks])) = 0;
  }


  /*------------------------------
    Check the files exist
    ------------------------------*/
  myprintf("Checking file existence\n");
  f1 = fopen(peakfilename,"rb");
  if (!f1) 
  {
    peaks_available=0; 
    myprintf("Could not open file [%s] for reading\n",peakfilename);
  }
  f2 = fopen(summitfilename,"rb");
  if (!f2) 
  {
    summits_available=0; 
    myprintf("Could not open file [%s] for reading\n",summitfilename);
  }
  if ((!f1) && (!f2)) 
    return(mappings);


  /*------------------------------
    Skip first lines (headers)
    ------------------------------*/
  myprintf("Skipping first %d peak lines\n",n_linestoskip_peaks);
  for (i=0;i<n_linestoskip_peaks;i++)
    readsuccess1=fgets(tline1,MAX_LINELENGTH,f1);

  myprintf("Skipping first %d summit lines\n",n_linestoskip_summits);
  for (i=0;i<n_linestoskip_summits;i++)
    readsuccess2=fgets(tline2,MAX_LINELENGTH,f1);


  /*------------------------------
    Read in the rest of the files
    ------------------------------*/
  chr_index=0;
  peakstarts=(int *)(mappings[chr_index][index_peakstart]);
  peakends=(int *)(mappings[chr_index][index_peakend]);
  peakscores=(float *)(mappings[chr_index][index_peakscore]);
  summitstarts=(int *)(mappings[chr_index][index_summitstart]);
  summitends=(int *)(mappings[chr_index][index_summitend]);
  summitheights=(float *)(mappings[chr_index][index_summitheight]);
  linenumbers=(int *)(mappings[chr_index][index_linenumber]);
  npeaks=*((int *)(mappings[chr_index][index_npeaks]));
  maxpeaks=*((int *)(mappings[chr_index][index_maxpeaks]));

  nlines=0;
  readsuccess1=fgets(tline1,MAX_LINELENGTH,f1);
  readsuccess2=fgets(tline2,MAX_LINELENGTH,f2);
  while ((readsuccess1 != NULL) || (readsuccess2 != NULL))
  {
    nlines=nlines+1;
    if ((nlines%100000)==0)
    {
      myprintf("Reading line %d\n",nlines);
      myprintf("%s", tline1);
    }

    if (readsuccess1)
    {

      /*------------------------------
	parse the peak line, 5 columns of interest
	------------------------------*/

      /* chromosome index 
      myprintf("Reading chromosome index\n"); */
    
      i1 = 0; 
      i2 = i1;
      while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 1 */
      memcpy(tempbuffer1,tline1+i1,i2-i1);
      tempbuffer1[i2-i1]='\0';
      /*myprintf("[%s]\n",tempbuffer1);*/
      line_chrindex1=-1;
      for (k=0;k<n_chromosomes;k++) 
	if (!strcmp(chromosomenames[k],tempbuffer1))
        {
	  line_chrindex1=k;
	  break;
	}
      if (line_chrindex1 < 0)
      {
	myprintf("Did not find a match for chromosome identifier [%s]\n", tempbuffer1);
	return(mappings);
      }

      /* myprintf("Read chromosome index %d\n", line_chrindex1);*/

      /* peak start 
      myprintf("Reading peak start\n");*/
    
      i1 = i2;
      while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 2 */
      i2 = i1;
      while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 2 */
      memcpy(tempbuffer1,tline1+i1,i2-i1);
      tempbuffer1[i2-i1]='\0';
      line_peakstart = atoi(tempbuffer1);
      /*myprintf("Read peak start %d\n", line_peakstart);*/

      /* peak end 
      myprintf("Reading peak end\n");*/
    
      i1 = i2;
      while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 3 */
      i2 = i1;
      while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 3 */
      memcpy(tempbuffer1,tline1+i1,i2-i1);
      tempbuffer1[i2-i1]='\0';
      line_peakend = atoi(tempbuffer1);
      /*myprintf("Read peak end %d\n", line_peakend);*/

      /* read name, we store this only for consistency checking
      myprintf("Reading peak name\n");*/
    
      i1 = i2;
      while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 4 */
      i2 = i1;
      while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 4 */
      memcpy(peakname1,tline1+i1,i2-i1);
      peakname1[i2-i1]='\0';

      /* peak score 
      myprintf("Reading peak score\n");*/
    
      i1 = i2;
      while ((tline1[i1] == ' ') || (tline1[i1] == 9)) i1=i1+1; /* start of column 5 */
      i2 = i1;
      while ((tline1[i2] != ' ') && (tline1[i2] != 9)) i2=i2+1; /* end of column 5 */
      memcpy(tempbuffer1,tline1+i1,i2-i1);
      tempbuffer1[i2-i1]='\0';
      line_peakscore = strtod(tempbuffer1, NULL);
      /*myprintf("Read peak score %f\n", line_peakscore);*/
    }

    if (readsuccess2)
    {
      /*------------------------------
	parse the summit line, 5 columns of interest
	------------------------------*/

      /* chromosome index 
      myprintf("Reading chromosome index\n"); */
    
      i1 = 0; 
      i2 = i1;
      while ((tline2[i2] != ' ') && (tline2[i2] != 9)) i2=i2+1; /* end of column 1 */
      memcpy(tempbuffer2,tline2+i1,i2-i1);
      tempbuffer2[i2-i1]='\0';
      /*myprintf("[%s]\n",tempbuffer2);*/
      line_chrindex2=-1;
      for (k=0;k<n_chromosomes;k++) 
	if (!strcmp(chromosomenames[k],tempbuffer2))
        {
	  line_chrindex2=k;
	  break;
	}
      if (line_chrindex2 < 0)
      {
	myprintf("Did not find a match for chromosome identifier [%s]\n", tempbuffer2);
	return(mappings);
      }

      /*myprintf("Read chromosome index %d\n", line_chrindex2);*/

      /* summit start 
      myprintf("Reading summit start\n");*/
    
      i1 = i2;
      while ((tline2[i1] == ' ') || (tline2[i1] == 9)) i1=i1+1; /* start of column 2 */
      i2 = i1;
      while ((tline2[i2] != ' ') && (tline2[i2] != 9)) i2=i2+1; /* end of column 2 */
      memcpy(tempbuffer2,tline2+i1,i2-i1);
      tempbuffer2[i2-i1]='\0';
      line_summitstart = atoi(tempbuffer2);
      /*myprintf("Read summit start %d\n", line_summitstart);*/

      /* summit end 
      myprintf("Reading summit end\n");*/
    
      i1 = i2;
      while ((tline2[i1] == ' ') || (tline2[i1] == 9)) i1=i1+1; /* start of column 3 */
      i2 = i1;
      while ((tline2[i2] != ' ') && (tline2[i2] != 9)) i2=i2+1; /* end of column 3 */
      memcpy(tempbuffer2,tline2+i1,i2-i1);
      tempbuffer2[i2-i1]='\0';
      line_summitend = atoi(tempbuffer2);
      /*myprintf("Read summit end %d\n", line_summitend);*/

      /* read name, we store this only for consistency checking
      myprintf("Reading summit name\n");*/
    
      i1 = i2;
      while ((tline2[i1] == ' ') || (tline2[i1] == 9)) i1=i1+1; /* start of column 4 */
      i2 = i1;
      while ((tline2[i2] != ' ') && (tline2[i2] != 9)) i2=i2+1; /* end of column 4 */
      memcpy(peakname2,tline2+i1,i2-i1);
      peakname2[i2-i1]='\0';

      /* summit height 
      myprintf("Reading summit height\n");*/
    
      i1 = i2;
      while ((tline2[i1] == ' ') || (tline2[i1] == 9)) i1=i1+1; /* start of column 5 */
      i2 = i1;
      while ((tline2[i2] != ' ') && (tline2[i2] != 9)) i2=i2+1; /* end of column 5 */
      memcpy(tempbuffer2,tline2+i1,i2-i1);
      tempbuffer2[i2-i1]='\0';
      line_summitheight = strtod(tempbuffer2, NULL);
      /*myprintf("Read summit height %f\n", line_summitheight);*/
    }

    /* Check consistency between peak and summit */
    if ((readsuccess1) && (readsuccess2))
    {
      if (line_chrindex1 != line_chrindex2)
      {
	myprintf("Mismatch in chromosome indices: peakfile %d vs summitfile %d\n", line_chrindex1, line_chrindex2);
	return(mappings);
      }
      if (strcmp(peakname1,peakname2))
      {
	myprintf("Mismatch in peak names: peakfile [%s] vs summitfile [%s]\n", peakname1, peakname2);
	return(mappings);
      }
    }

    if ((readsuccess1) && (!readsuccess2))
      line_chrindex2 = line_chrindex1;
    if ((readsuccess2) && (!readsuccess1))
      line_chrindex1 = line_chrindex2;


    /* start working on the chromosome identified on the line, if we weren't already working on it */
    /*myprintf("Starting to work on correct chromosome\n");*/
    
    if (line_chrindex1 != chr_index)
    {
      /*myprintf("before: chr %d, npeaks %d, maxpeaks %d",chr_index,npeaks,maxpeaks);*/
      mappings[chr_index][index_peakstart] = peakstarts;
      mappings[chr_index][index_peakend] = peakends;
      mappings[chr_index][index_peakscore] = peakscores;
      mappings[chr_index][index_summitstart] = summitstarts;
      mappings[chr_index][index_summitend] = summitends;
      mappings[chr_index][index_summitheight] = summitheights;
      mappings[chr_index][index_linenumber] = linenumbers;
      *((int *)(mappings[chr_index][index_npeaks])) = npeaks;
      *((int *)(mappings[chr_index][index_maxpeaks])) = maxpeaks;
    
      chr_index = line_chrindex1;
      peakstarts=(int *)(mappings[chr_index][index_peakstart]);
      peakends=(int *)(mappings[chr_index][index_peakend]);
      peakscores=(float *)(mappings[chr_index][index_peakscore]);
      summitstarts=(int *)(mappings[chr_index][index_summitstart]);
      summitends=(int *)(mappings[chr_index][index_summitend]);
      summitheights=(float *)(mappings[chr_index][index_summitheight]);
      linenumbers=(int *)(mappings[chr_index][index_linenumber]);
      npeaks=*((int *)(mappings[chr_index][index_npeaks]));
      maxpeaks=*((int *)(mappings[chr_index][index_maxpeaks]));
      /*myprintf("after: chr %d, npeaks %d, maxpeaks %d",chr_index,npeaks,maxpeaks);*/
    }

    /* add the peak/summit to the data of the chromosome */
    /*myprintf("Adding the lines to the data for the chromosome %d, %d peaks, %d maxpeaks\n", chr_index, npeaks,maxpeaks);*/
    
    npeaks=npeaks+1;    
    if (npeaks>maxpeaks)
    {
      maxpeaks=npeaks*2;
      peakstarts=(int *)myrealloc(peakstarts, sizeof(int)*maxpeaks);
      peakends=(int *)myrealloc(peakends, sizeof(int)*maxpeaks);
      peakscores=(float *)myrealloc(peakscores, sizeof(float)*maxpeaks);
      summitstarts=(int *)myrealloc(summitstarts, sizeof(int)*maxpeaks);
      summitends=(int *)myrealloc(summitends, sizeof(int)*maxpeaks);
      summitheights=(float *)myrealloc(summitheights, sizeof(float)*maxpeaks);
      linenumbers=(int *)myrealloc(linenumbers, sizeof(int)*maxpeaks);
    }
  
    if (readsuccess1)
    {
      peakstarts[npeaks-1]=line_peakstart;
      peakends[npeaks-1]=line_peakend;
      peakscores[npeaks-1]=line_peakscore;
    }
    else
    {
      peakstarts[npeaks-1]=-1;
      peakends[npeaks-1]=-1;
      peakscores[npeaks-1]=-1;
    }

    if (readsuccess2)
    {
      summitstarts[npeaks-1]=line_summitstart;
      summitends[npeaks-1]=line_summitend;
      summitheights[npeaks-1]=line_summitheight;
    }
    else
    {
      summitstarts[npeaks-1]=-1;
      summitends[npeaks-1]=-1;
      summitheights[npeaks-1]=-1;
    }

    linenumbers[npeaks-1]=nlines;  /* line number where this peak occurred in the files */

    /* read in the next lines */
    /*myprintf("Reading next lines\n");*/
    readsuccess1=fgets(tline1,MAX_LINELENGTH,f1);
    readsuccess2=fgets(tline2,MAX_LINELENGTH,f2);

    /* debug early stopping point */
    /* if (nlines > 1000) readsuccess = NULL; */
  }

  myprintf("Finished reading file, storing working data\n");
  

  /* store the data of the last chromosome we were working on */
  mappings[chr_index][index_peakstart] = peakstarts;
  mappings[chr_index][index_peakend] = peakends;
  mappings[chr_index][index_peakscore] = peakscores;
  mappings[chr_index][index_summitstart] = summitstarts;
  mappings[chr_index][index_summitend] = summitends;
  mappings[chr_index][index_summitheight] = summitheights;
  mappings[chr_index][index_linenumber] = linenumbers;
  *((int *)(mappings[chr_index][index_npeaks])) = npeaks;
  *((int *)(mappings[chr_index][index_maxpeaks])) = maxpeaks;

  myprintf("Discarding empty buffer data\n");


  /* for all chromosomes, discard empty data */
  for (chr_index=0; chr_index<n_chromosomes;chr_index++)
  {
    peakstarts=(int *)(mappings[chr_index][index_peakstart]);
    peakends=(int *)(mappings[chr_index][index_peakend]);
    peakscores=(float *)(mappings[chr_index][index_peakscore]);
    summitstarts=(int *)(mappings[chr_index][index_summitstart]);
    summitends=(int *)(mappings[chr_index][index_summitend]);
    summitheights=(float *)(mappings[chr_index][index_summitheight]);
    linenumbers=(int *)(mappings[chr_index][index_linenumber]);
    npeaks=*((int *)(mappings[chr_index][index_npeaks]));
    maxpeaks=*((int *)(mappings[chr_index][index_maxpeaks]));

    maxpeaks=npeaks;
    peakstarts=(int *)myrealloc(peakstarts, sizeof(int)*maxpeaks);
    peakends=(int *)myrealloc(peakends, sizeof(int)*maxpeaks);
    peakscores=(float *)myrealloc(peakscores, sizeof(float)*maxpeaks);
    summitstarts=(int *)myrealloc(summitstarts, sizeof(int)*maxpeaks);
    summitends=(int *)myrealloc(summitends, sizeof(int)*maxpeaks);
    summitheights=(float *)myrealloc(summitheights, sizeof(float)*maxpeaks);
    linenumbers=(int *)myrealloc(linenumbers, sizeof(int)*maxpeaks);

    mappings[chr_index][index_peakstart] = peakstarts;
    mappings[chr_index][index_peakend] = peakends;
    mappings[chr_index][index_peakscore] = peakscores;
    mappings[chr_index][index_summitstart] = summitstarts;
    mappings[chr_index][index_summitend] = summitends;
    mappings[chr_index][index_summitheight] = summitheights;
    mappings[chr_index][index_linenumber] = linenumbers;
    *((int *)(mappings[chr_index][index_npeaks])) = npeaks;
    *((int *)(mappings[chr_index][index_maxpeaks])) = maxpeaks;
  }

  myprintf("Closing files\n");

  fclose(f1);
  fclose(f2);

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
  mxArray *tempfilename1;
  mxArray *tempfilename2;
  char filename1[MAX_LINELENGTH];
  char filename2[MAX_LINELENGTH];

  mxArray *tempmappings;
  mxArray *temppeakstarts;
  mxArray *temppeakends;
  mxArray *temppeakscores;
  mxArray *tempsummitstarts;
  mxArray *tempsummitends;
  mxArray *tempsummitheights;
  mxArray *templinenumbers;
  mxArray *tempnpeaks;
  mxArray *tempmaxpeaks;
  int *temppeakstarts_data;
  int *temppeakends_data;
  float *temppeakscores_data;
  int *tempsummitstarts_data;
  int *tempsummitends_data;
  float *tempsummitheights_data;
  int *templinenumbers_data;
  int *tempnpeaks_data;
  int *tempmaxpeaks_data;

  void ***mappings;
  int npeaks;
  int maxpeaks;
  int *peakstarts;
  int *peakends;
  float *peakscores;
  int *summitstarts;
  int *summitends;
  float *summitheights;
  int *linenumbers;

  int chr_index, j, k;
  int tempindex[2];

  tempfilename1 = prhs[0];
  mxGetString(tempfilename1, filename1, MAX_LINELENGTH);

  tempfilename2 = prhs[1];
  mxGetString(tempfilename2, filename2, MAX_LINELENGTH);

  myprintf("Trying to read from peak file [%s] and summit file [%s] \n",filename1,filename2);
  mappings = read_peak_summit_files(filename1,filename2);

  tempmappings = mxCreateCellMatrix(n_chromosomes,9);

  for (chr_index=0;chr_index<n_chromosomes;chr_index++)
  {    
    peakstarts=(int *)(mappings[chr_index][index_peakstart]);
    peakends=(int *)(mappings[chr_index][index_peakend]);
    peakscores=(float *)(mappings[chr_index][index_peakscore]);
    summitstarts=(int *)(mappings[chr_index][index_summitstart]);
    summitends=(int *)(mappings[chr_index][index_summitend]);
    summitheights=(float *)(mappings[chr_index][index_summitheight]);
    linenumbers=(int *)(mappings[chr_index][index_linenumber]);
    npeaks=*((int *)(mappings[chr_index][index_npeaks]));
    maxpeaks=*((int *)(mappings[chr_index][index_maxpeaks]));

    temppeakstarts = mxCreateNumericMatrix(npeaks,1,mxINT32_CLASS,0);
    temppeakstarts_data = (int *)mxGetData(temppeakstarts);

    temppeakends = mxCreateNumericMatrix(npeaks,1,mxINT32_CLASS,0);
    temppeakends_data = (int *)mxGetData(temppeakends);

    temppeakscores = mxCreateNumericMatrix(npeaks,1,mxSINGLE_CLASS,0);
    temppeakscores_data = (float *)mxGetData(temppeakscores);

    tempsummitstarts = mxCreateNumericMatrix(npeaks,1,mxINT32_CLASS,0);
    tempsummitstarts_data = (int *)mxGetData(tempsummitstarts);

    tempsummitends = mxCreateNumericMatrix(npeaks,1,mxINT32_CLASS,0);
    tempsummitends_data = (int *)mxGetData(tempsummitends);

    tempsummitheights = mxCreateNumericMatrix(npeaks,1,mxSINGLE_CLASS,0);
    tempsummitheights_data = (float *)mxGetData(tempsummitheights);

    templinenumbers = mxCreateNumericMatrix(npeaks,1,mxINT32_CLASS,0);
    templinenumbers_data = (int *)mxGetData(templinenumbers);

    tempnpeaks = mxCreateNumericMatrix(1,1,mxINT32_CLASS,0);
    tempnpeaks_data = (int *)mxGetData(tempnpeaks);

    tempmaxpeaks = mxCreateNumericMatrix(1,1,mxINT32_CLASS,0);
    tempmaxpeaks_data = (int *)mxGetData(tempmaxpeaks);

    tempnpeaks_data[0] = npeaks;
    tempmaxpeaks_data[0] = maxpeaks;
    for (j=npeaks-1;j>=0;j--)
    {
      temppeakstarts_data[j] = peakstarts[j];
      temppeakends_data[j] = peakends[j];
      temppeakscores_data[j] = peakscores[j];
      tempsummitstarts_data[j] = summitstarts[j];
      tempsummitends_data[j] = summitends[j];
      tempsummitheights_data[j] = summitheights[j];
      templinenumbers_data[j] = linenumbers[j];
    }

    tempindex[0]=chr_index;
    tempindex[1]=index_peakstart;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,temppeakstarts);

    tempindex[0]=chr_index;
    tempindex[1]=index_peakend;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,temppeakends);

    tempindex[0]=chr_index;
    tempindex[1]=index_peakscore;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,temppeakscores);

    tempindex[0]=chr_index;
    tempindex[1]=index_summitstart;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempsummitstarts);

    tempindex[0]=chr_index;
    tempindex[1]=index_summitend;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempsummitends);

    tempindex[0]=chr_index;
    tempindex[1]=index_summitheight;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempsummitheights);

    tempindex[0]=chr_index;
    tempindex[1]=index_linenumber;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,templinenumbers);

    tempindex[0]=chr_index;
    tempindex[1]=index_npeaks;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempnpeaks);

    tempindex[0]=chr_index;
    tempindex[1]=index_maxpeaks;
    j = mxCalcSingleSubscript(tempmappings,2,tempindex);
    mxSetCell(tempmappings,j,tempmaxpeaks);


    mxFree(peakstarts);
    mxFree(peakends);
    mxFree(peakscores);
    mxFree(summitstarts);
    mxFree(summitends);
    mxFree(summitheights);
    mxFree(linenumbers);
  }

  plhs[0] = tempmappings;
  
}


