beep_on_error(1);
warning("error",'Octave:divide-by-zero');


% mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
% h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
% pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/'
% analysisdir='/Users/hasanogul/jaakkos_files/synergy/analyses/'

mybasedir='~/synergy_data/tempcodebranch/';
h3k4me3dir='~/synergy_data/H3K4me3/Mapping_results/'
pol2dir='~/synergy_data/PolII/Mapping_results/'
rnadir='~/synergy_data/RNA/Mapping_results/'
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




cd(pol2dir)
load pol2_for_matti_ver3.mat   % provides variable bininfo

load all_gene_pol2bins.mat  % provides: pol2bins 
load all_empty_pol2bins.mat % provides: allemptybins emptybininfo
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
% Read in exon information
%-----------------------------------------------------------------

cd(rnadir);
exonfile=read_stringfile('exons.txt',[32 9 10 13]);
firstexonstarts=ones(size(bininfo,1),1)*nan;
for k=1:(length(exonfile)-1)/4,
  ensemblid=str2double(exonfile{(k-1)*4+2}{1}(5:end));
  for l=1:length(exonfile{(k-1)*4+3}),
    exonstarts(l)=str2double(exonfile{(k-1)*4+4}{l});
  end;
  for l=1:length(exonfile{(k-1)*4+4}),
    exonends(l)=str2double(exonfile{(k-1)*4+5}{l});
  end;
  % all we actually need for now is the start of the first exon
  firstexonstart=min(exonstarts);
  
  % find the index in bininfo that corresponds to this ensembl-id
  gene_index=find(bininfo(:,5)==ensemblid);
  if length(gene_index)==1,
    firstexonstarts(gene_index)=firstexonstart;
  end;     
end;

% find the bin corresponding to the exon start
firstexonbins=ones(size(bininfo,1),1)*nan;
nextexonbins=ones(size(bininfo,1),1)*nan;
for k=1:size(bininfo,1),
  if ~isnan(firstexonstarts(k)),
    % start of the first exon -30bp
    differencefromend=bininfo(k,3)-(firstexonstarts(k)-30);
    firstexonbins(k)=1+floor(differencefromend/200);
    if firstexonbins(k)>binlengths(k), firstexonbins(k)=binlengths(k);end;
    
    % that location + 330bp
    differencefromend=bininfo(k,3)-(firstexonstarts(k)+300);
    if differencefromend < -200,
      nextexonbins(k)=-100;
    elseif differencefromend < 0,
      nextexonbins(k)=0;
    else      
      nextexonbins(k)=1+floor(differencefromend/200);
    end;
    if nextexonbins(k)==0, nextexonbins(k)=1;end;
    
  else
    % just take the start of the gene
    firstexonbins(k)=binlengths(k);
    if firstexonbins(k)>binlengths(k), firstexonbins(k)=binlengths(k);end;

    % that location + 300bp
    differencefromend=bininfo(k,3)-(bininfo(k,2)+300);
    if differencefromend < -200,
      nextexonbins(k)=-100;
    elseif differencefromend < 0,
      nextexonbins(k)=0;
    else      
      nextexonbins(k)=1+floor(differencefromend/200);
    end;
    if nextexonbins(k)==0, nextexonbins(k)=1;end;
  end;  
end;




%-----------------------------------------------------------------
% Perform noise thresholding
%-----------------------------------------------------------------

% For each time point, compute a noise level (average activity per
% bin in hand-picked 'empty regions'. Note that the empty regions
% have been hand-picked based on POL2 activity from chromosome 1 only.
emptysums=zeros(n_timepoints,1);nemptybins=zeros(n_timepoints,1);
for k=1:n_timepoints,
  k
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
    l
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
  if mod(k,1000)==0,
    k
  end;
  
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
  if mod(i,1000)==0,
    i
  end;
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
% polsummaryoffsetfraction=0;
pol2_summaryseries_start=nan*ones(size(pol2bins,1),n_timepoints);
pol2_summaryseries_middle=nan*ones(size(pol2bins,1),n_timepoints);
pol2_summaryseries_end=nan*ones(size(pol2bins,1),n_timepoints);
for i=1:size(pol2bins,1),
  if mod(i,100)==0,
    i
  end;
  
  % summarize by first bins 
  % (note that first bins = bins closest to transcription end)
  polsummarylength=ceil(binlengths(i)*polsummaryfraction);
  % polsummaryoffset=ceil(binlengths(i)*polsummaryoffsetfraction);

  start_begin=firstexonbins(i);
  middle_begin=nextexonbins(i);
  end_begin=polsummarylength;
  start_length=start_begin-middle_begin;
  middle_length=middle_begin-end_begin;
  end_length=end_begin;
  
  if (binlengths(i)>=start_begin) && (start_begin>middle_begin) && (middle_begin>end_begin) && (end_begin>=1),
    for j=1:n_timepoints,
      pol2_summaryseries_start(i,j)=sum(pol2bins{i,j}(middle_begin+1:start_begin))/(start_length*200);
      pol2_summaryseries_middle(i,j)=sum(pol2bins{i,j}(end_begin+1:middle_begin))/(middle_length*200);
      pol2_summaryseries_end(i,j)=sum(pol2bins{i,j}(1:end_begin))/(end_length*200);
    end;
  end;  
end;


%---------------------------------------------------
% normalize H3K4me3 time points by noise levels estimated
% from geommean procedure earlier.
%---------------------------------------------------
for j=1:n_timepoints,
  pol2_summaryseries_start(:,j)=pol2_summaryseries_start(:,j)/tempmedians(j);
  pol2_summaryseries_middle(:,j)=pol2_summaryseries_middle(:,j)/tempmedians(j);
  pol2_summaryseries_end(:,j)=pol2_summaryseries_end(:,j)/tempmedians(j);
end;
                    

beep();
