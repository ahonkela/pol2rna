beep_on_error(1);
warning("error",'Octave:divide-by-zero');


%mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
%h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
%pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/'
%analysisdir='/Users/hasanogul/jaakkos_files/synergy/analyses/'

mybasedir='~/synergy_data/tempcodebranch/';
h3k4me3dir='~/synergy_data/H3K4me3/Mapping_results/'
pol2dir='~/synergy_data/PolII/Mapping_results/'
groseqdir='~/synergy_data/GROseq/Mapping_results/'
analysisdir='~/synergy_data/analyses/'


% for kernel-level computations
path1=[mybasedir 'kern/matlab/jaakko_testversion/']
% for model-level computations
path2=[mybasedir 'gpsim/matlab/jaakko_testversion/'];
% for optimiDefaultConstraint.m
path3=[mybasedir 'optimi/matlab/jaakko_testversion/'];
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




cd(groseqdir)
load all_gene_groseqbins.mat  % provides: groseqbins bininfo
% load all_empty_groseqbins.mat % provides: allemptybins emptybininfo
n_timepoints=8;

%  bininfo(i,1)=line_chrindex;                    % chromosome index
%  bininfo(i,2)=str2double(genetemp{3}(2:end-1)); % start location
%  bininfo(i,3)=str2double(genetemp{4}(2:end-1)); % end location
%  bininfo(i,4)=i;                                % file line number
%  bininfo(i,5)=geneids(i);                       % ensembl id
%  bininfo(i,6)=genestrands(i);                   % strand sign (+1 or -1)


binlengths=zeros(size(groseqbins,1),1);
for k=1:size(groseqbins,1),
  binlengths(k)=length(groseqbins{k,1});
end;



%-----------------------------------------------------------------
% Perform noise thresholding
%-----------------------------------------------------------------

apply_noisethreshold=0;
if apply_noisethreshold,
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

  %---------------------------------------------------
  % Substract noise level, 
  % thresholding values below noise level to zero. 
  %---------------------------------------------------
  for l=1:n_timepoints,
    for k=1:size(groseqbins,1),
      groseqbins{k,l}=groseqbins{k,l}-noiselevels(l);
      I=find(groseqbins{k,l}<0);
      groseqbins{k,l}(I)=0;
    end;
  end;
end;


%-----------------------------------------------------------------
% Compute geometric mean - median based normalization factors for
% time points
%-----------------------------------------------------------------

% find genes that have a sufficient level above noise in all time points
genemeans=zeros(size(groseqbins,1),n_timepoints);
for k=1:size(groseqbins,1),
  for l=1:n_timepoints,
    if length(groseqbins{k,l})>0,
      genemeans(k,l)=mean(groseqbins{k,l});
    end;
  end;
end;
readthreshold=(200/200)*200; % GROseq produces much lower readcounts than
                     % ChIPseq, so needs a lower threshold. This
                     % one is just arbitrarily handpicked: there
                     % must be at least five single-basepair long
                     % reads among the 200 basepairs in each bin.
enoughdata=genemeans>readthreshold;

% compute geommeanmedian normalization from the above selected genes, using sum over whole gene
tempgeommeans=zeros(size(groseqbins,1),1);
for i=1:size(groseqbins,1),
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
% Create GROseq summary series
%-----------------------------------------------------------------

% Version 3: sum of last few bins (fixed fraction of gene length, 
% towards the end of the gene), normalized by total length of those bins                    

%groseqsummaryfraction=0.2;
%groseqsummaryoffsetfraction=0;
groseq_summaryseries=nan*ones(size(groseqbins,1),n_timepoints);
for i=1:size(groseqbins,1),
  % summarize by last bins 
  % (note that last bins = bins closest to transcription start)
  groseqsummaryoffset=0;
  
  % use the same summarization as in the Hah et al. Cell paper
  % "A rapid, extensive, and transient transcriptional response to
  % estrogen signaling in breast cancer cells". Note that bins are
  % 200 basepairs long (except for possibly the last bin which can
  % be shorter).
  
  % these are zero-based indices from the gene start
  startposition_from_genestart = 5;
  endposition_from_genestart = 13*5 -1;
  if (endposition_from_genestart > binlengths(i)-1),
    endposition_from_genestart = binlengths(i)-1;
  end;
  
  % translate the positions to positions from gene end...
  groseqsummaryoffset = binlengths(i)-1-endposition_from_genestart;
  groseqsummarylength = endposition_from_genestart-startposition_from_genestart+1;
  
  if (groseqsummarylength>0),
    % make sure we are not taking the peak at the transcription start...
    if (binlengths(i)>=groseqsummarylength+groseqsummaryoffset),
      for j=1:n_timepoints,
        groseq_summaryseries(i,j)=...
            sum(groseqbins{i,j}((1+groseqsummaryoffset):(groseqsummarylength+groseqsummaryoffset))) ...
            / groseqsummarylength*200;   % normalize by total length of the chosen bins
      end;
    end; 
  end;
end;



%---------------------------------------------------
% normalize GROseq time points by noise levels estimated
% from geommean procedure earlier.
%---------------------------------------------------
for j=1:n_timepoints,
  groseq_summaryseries(:,j)=groseq_summaryseries(:,j)/tempmedians(j);
end;
                    

beep();
