%beep_on_error(1);
%warning('error','Octave:divide-by-zero');


mybasedir='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
% h3k4me3dir='/share/mi/workspace/jtpelto/synergy/synergy_data/H3K4me3/Mapping_results/'
pol2dir='/share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed/'
% analysisdir='/share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed/'

% for kernel-level computations
path1=[mybasedir 'kern/matlab/']
% for model-level computations
path2=[mybasedir 'gpsim/matlab/'];
% for optimiDefaultConstraint.m
path3=[mybasedir 'optimi/matlab/'];
% for lnDiffErfs.m
path4=[mybasedir 'ndlutil/matlab/'];
% for addPrior.m
path5=[mybasedir 'prior/matlab/'];
% for dist2.m
path6=[mybasedir 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m
path7=[mybasedir 'mltools/matlab/'];
% for various experiment things
path8=[mybasedir 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)




cd(pol2dir)
load all_gene_pol2bins_2012_03.mat  % provides: pol2bins bininfo
load all_empty_pol2bins_2012_03.mat % provides: allemptybins emptybininfo
n_timepoints=10;

%  bininfo(i,1)=line_chrindex;                    % chromosome index
%  bininfo(i,2)=str2double(genetemp{3}(2:end-1)); % start location
%  bininfo(i,3)=str2double(genetemp{4}(2:end-1)); % end location
%  bininfo(i,4)=i;                                % file line number
%  bininfo(i,5)=geneids(i);                       % ensembl id
%  bininfo(i,6)=genestrands(i);                   % strand sign (+1 or -1)


binlengths=zeros(size(pol2bins,1),1);
for k=1:size(pol2bins,1),
  binlengths(k)=length(pol2bins{k,1});
end;



%-----------------------------------------------------------------
% Perform noise thresholding
%-----------------------------------------------------------------

% For each time point, compute a noise level (average activity per
% bin in hand-picked 'empty regions'. Note that the empty regions
% have been hand-picked based on POL2 activity from chromosome 1 only.
emptysums=zeros(n_timepoints,1);nemptybins=zeros(n_timepoints,1);
for k=1:n_timepoints,
  for l=1:74,
    emptysums(k)=emptysums(k)+sum(allemptybins{l,k});
    nemptybins(k)=nemptybins(k)+length(allemptybins{l,k});
  end;
end;
noiselevels=emptysums./nemptybins;


apply_noisethreshold=1;
if apply_noisethreshold,
  %---------------------------------------------------
  % Substract noise level, 
  % thresholding values below noise level to zero. 
  %---------------------------------------------------
  for l=1:n_timepoints,
    for k=1:size(pol2bins,1),
      pol2bins{k,l}=pol2bins{k,l}-noiselevels(l);
      I=find(pol2bins{k,l}<0);
      pol2bins{k,l}(I)=0;
    end;
  end;
end;


%-----------------------------------------------------------------
% Compute geometric mean - median based normalization factors for
% time points
%-----------------------------------------------------------------

% find genes that have a sufficient level above noise in all time points
genemeans=zeros(size(pol2bins,1),n_timepoints);
for k=1:size(pol2bins,1),
  for l=1:n_timepoints,
    if length(pol2bins{k,l})>0,
      genemeans(k,l)=mean(pol2bins{k,l});
    end;
  end;
end;
readthreshold=5*200;
enoughdata=genemeans>readthreshold;

% compute geommeanmedian normalization from the above selected genes, using sum over whole gene
tempgeommeans=zeros(size(pol2bins,1),1);
for i=1:size(pol2bins,1),
  % find the time points of this gene that have enough data,
  % excluding the first time point
  Igeommean=find(enoughdata(i,2:n_timepoints)==1)+1;
  
  % compute geometric mean of the selected time points for this gene
  if length(Igeommean>0),
    tempgeommeans(i)=exp(mean(log(genemeans(i,Igeommean))));  
  end;
end;
% compute the normalization factors for the time points as a median
% normalization factor across genes.
Imedians=find(tempgeommeans>0);
tempmedians=zeros(1,n_timepoints);
tempmedians(1)=1;
for j=2:n_timepoints,
  tempmedians(j)=median(genemeans(Imedians,j)./tempgeommeans(Imedians));
end;



%-----------------------------------------------------------------
% Create POL2 summary series
%-----------------------------------------------------------------

% Version 3: sum of last few bins (fixed fraction of gene length, 
% towards the end of the gene), normalized by total length of those bins                    
polsummaryfraction=0.2;
polsummaryoffsetfraction=0;
pol2_summaryseries=nan*ones(size(pol2bins,1),n_timepoints);
pol2_ndata=nan*ones(size(pol2bins,1),n_timepoints);
pol2_sums=nan*ones(size(pol2bins,1),n_timepoints);
normalizers=nan*ones(size(pol2bins,1),n_timepoints);
for i=1:size(pol2bins,1),
  % summarize by first bins 
  % (note that first bins = bins closest to transcription end)
  polsummarylength=ceil(binlengths(i)*polsummaryfraction);
  polsummaryoffset=ceil(binlengths(i)*polsummaryoffsetfraction);

  if (polsummarylength>0),
    % make sure we are not taking the peak at the transcription start...
    if (binlengths(i)>=polsummarylength+polsummaryoffset),
      for j=1:n_timepoints,
        pol2_summaryseries(i,j)=...
            sum(pol2bins{i,j}((1+polsummaryoffset):(polsummarylength+polsummaryoffset))) ...
            / (polsummarylength*200);   % normalize by total length of the chosen bins
        pol2_ndata(i,j)=...
            sum(pol2bins{i,j}((1+polsummaryoffset):(polsummarylength+polsummaryoffset))) ...
            / (polsummarylength*200);   % normalize by total length of the chosen bins
        pol2_sums(i,j)=...
            sum(pol2bins{i,j}((1+polsummaryoffset):(polsummarylength+polsummaryoffset)));
	normalizers(i,j)=polsummarylength*200;
      end;
    end; 
  end;
end;


%---------------------------------------------------
% normalize H3K4me3 time points by noise levels estimated
% from geommean procedure earlier.
%---------------------------------------------------
for j=1:n_timepoints,
  pol2_summaryseries(:,j)=pol2_summaryseries(:,j)/tempmedians(j);
end;
                    

beep();

% save pol2_summaryseries_2012_03.mat pol2_summaryseries -mat