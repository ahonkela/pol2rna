% mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
% h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
% pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/'
% analysisdir='/Users/hasanogul/jaakkos_files/synergy/analyses/'

mybasedir='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
pol2gzdir = '/share/synergy/data/2012-03_PolII/Mapping_results/';
pol2dir='/share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed/'
analysisdir='/share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed/'




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

  
  
  
  
  
  cd(pol2dir)
  
  
  %-----------------------------
  % For each time point, count amount of POL2-read-basepairs that 
  % overlap each bin, weighted by scores of the reads.
  %-----------------------------
  
  %allbins=zeros(n_genes,10);
  filenames={'pol0_2012_03.mat','pol5_2012_03.mat','pol10_2012_03.mat','pol20_2012_03.mat','pol40_2012_03.mat','pol80_2012_03.mat','pol160_2012_03.mat','pol320_2012_03.mat','pol640_2012_03.mat','pol1280_2012_03.mat'};

  % dvalues=[230 214 219 213 217 223 215 208 213 212]; % dvalues obtained from MACS output
  dvalues=[0 0 0 0 0 0 0 0 0 0]; % do not shift or extend the reads
  
  ntimepoints=length(filenames);
  % load up all data
  allgenebins=cell(size(emptybininfo,1),ntimepoints);
  fprintf(1,'Initializing allgenebins\n');
%  load all_gene_h3k4me3bins.mat
%  allgenebins=h3k4me3bins;
  
  for timepoint=1:ntimepoints,
  %for timepoint=9:ntimepoints,
    fprintf(1,'Loading data of time point %d\n',timepoint);
    load(filenames{timepoint});  % provides variable temppol
    d=dvalues(timepoint);

    % For the 2012-03 data, MACS seems to keep only 1 duplicate;
    % let's do the same here even if we do not use MACS for shifting/lengthening reads
    max_duplicates=1;  

    subbin_length=200;
    %binsthistime=zeros(n_genes,1);
    
    for chr_index=1:n_chromosomes,
      timepoint
      chr_index

      % ensure that read scores are in double format, not single      
      if chr_index<=size(temppol,1),
        temppol{chr_index,4}=double(temppol{chr_index,4});
      end;
      
      I=find(emptybininfo(:,1)==chr_index);
      tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(emptybininfo(I,2)),int32(emptybininfo(I,3)),int8(emptybininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allgenebins{I(k),timepoint}=tempbins{k};
      end;
    end;
    clear temppol;
  end;
  allemptybins = allgenebins;
  save all_empty_pol2bins_2012_09.mat allemptybins emptybininfo -mat
  
end;
