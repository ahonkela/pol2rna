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

mexcodedir=path8;



if 1,
  cd(mexcodedir)
  mex compute_pol2activityovergenes_c.c

  chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
  n_chromosomes=length(chromosomenames);
  
  
  %-----------------------------
  % Read in empty region information
  %-----------------------------

  % read empty region chromosome ids and start and end locations
  emptylocs=read_stringfile([analysisdir 'empty_regions.txt'],[32 9 10 13]);
  n_emptyregions=length(emptylocs);


  %-----------------------------
  % Create bin array, one bin for each empty region
  %-----------------------------
  emptybininfo=zeros(length(emptylocs),6);
  for i=1:length(emptylocs),
    emptytemp=emptylocs{i};
    
    % get chromosome index in chromosome list
    emptychromosome=emptytemp{1}(1:end);
    line_chrindex=-1;
    for k=1:n_chromosomes,
      if strcmp(emptychromosome,chromosomenames{k})==1,
        line_chrindex=k;
        break;
      end;
    end;
    
    emptybininfo(i,1)=line_chrindex; % chromosome index
    emptybininfo(i,2)=str2double(emptytemp{2}(1:end)); % start location
    emptybininfo(i,3)=str2double(emptytemp{3}(1:end)); % end location
    emptybininfo(i,4)=i;
    emptybininfo(i,5)=i;  % just a dummy ID for empty regions
    emptybininfo(i,6)=1;  % just a dummy strand ID for empty regions
  end;



  %-----------------------------
  % Sort bin array in ascending order of start location. 
  % The POL2 activity computation code assumes bins have been sorted like that.
  %-----------------------------
  [y,I1]=sort(emptybininfo(:,2));
  emptybininfo=emptybininfo(I1,:);

  
  
  
  
  
  cd(groseqdir)
  
  
  %-----------------------------
  % For each time point, count amount of POL2-read-basepairs that 
  % overlap each bin, weighted by scores of the reads.
  %-----------------------------
  
  %allbins=zeros(n_genes,10);
  filenames={'groseq_0a.mat','groseq_0b.mat','groseq_10a.mat','groseq_10b.mat','groseq_40a.mat','groseq_40b.mat','groseq_160a.mat','groseq_160b.mat'};
  dvalues=[0 0 0 0 0 0 0 0]; % dvalues obtained from MACS output

  
  ntimepoints=length(filenames);
  % load up all data
  allgenebins=cell(size(emptybininfo,1),ntimepoints);
  fprintf(1,'Initializing allgenebins\n');
  
  for timepoint=1:ntimepoints,
  %for timepoint=9:ntimepoints,
    fprintf(1,'Loading data of time point %d\n',timepoint);
    load(filenames{timepoint});  % provides variable tempgroseq
    d=dvalues(timepoint);
    max_duplicates=2;  
    subbin_length=200;
    %binsthistime=zeros(n_genes,1);

    % set read scores to uniform 1 for groseq data
    for k=1:size(tempgroseq,1),
      tempgroseq{k,4}(:)=1;
    end;    
    
    for chr_index=1:n_chromosomes,
      timepoint
      chr_index

      % ensure that read scores are in double format, not single      
      if chr_index<=size(tempgroseq,1),
        tempgroseq{chr_index,4}=double(tempgroseq{chr_index,4});
      end;
      
      I=find(emptybininfo(:,1)==chr_index);
      tempbins=compute_pol2activityovergenes_c(tempgroseq,int32(chr_index),int32(length(I)),int32(emptybininfo(I,2)),int32(emptybininfo(I,3)),int8(emptybininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allgenebins{I(k),timepoint}=tempbins{k};
      end;
    end;
    clear tempgroseq;
  end;
  allemptybins = allgenebins;
  save all_empty_groseqbins.mat allemptybins emptybininfo -mat
  
end;
