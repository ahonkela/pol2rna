mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/'
analysisdir='/Users/hasanogul/jaakkos_files/synergy/analyses/'

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



cd(mexcodedir)
mex compute_pol2activityovergenes_c.c
cd(h3k4me3dir)

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
n_chromosomes=length(chromosomenames);

n_timepoints=9;

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



%-----------------------------
% For each time point, count amount of POL2-read-basepairs that 
% overlap each bin, weighted by scores of the reads.
%-----------------------------

filenames={'h3k4me3_0.mat','h3k4me3_10.mat','h3k4me3_20.mat','h3k4me3_40.mat','h3k4me3_80.mat','h3k4me3_160.mat','h3k4me3_320.mat','h3k4me3_640.mat','h3k4me3_1280.mat'};
dvalues=[192 192 196 196 180 197 186 198 211]; % dvalues obtained from MACS output


% load up all data
allemptybins=cell(size(emptybininfo,1),n_timepoints);
for timepoint=1:9,
  load(filenames{timepoint});  % this provides the variable
                               % temph3k4me3, containing bins for
                               % that time point.
  d=dvalues(timepoint);
  max_duplicates=2;  
  subbin_length=200;
  
  for chr_index=1:n_chromosomes,
    timepoint
    chr_index
    I=find(emptybininfo(:,1)==chr_index);
    if length(I)>0,
      tempbins=compute_pol2activityovergenes_c(temph3k4me3,int32(chr_index),int32(length(I)),int32(emptybininfo(I,2)),int32(emptybininfo(I,3)),int8(emptybininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allemptybins{I(k),timepoint}=tempbins{k};
      end;
    end;
  end;
  clear temppol;
end;

save all_empty_h3k4me3bins.mat allemptybins emptybininfo -mat


