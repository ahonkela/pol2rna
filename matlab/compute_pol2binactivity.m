addpath ~/mlprojects/pol2rnaseq/matlab
cd ~/mlprojects/pol2rnaseq/matlab/
mex compute_pol2activityovergenes_c.c
cd ~/synergy_data/PolII/Mapping_results

chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
n_chromosomes=length(chromosomenames);


%-----------------------------
% Read in gene information
%-----------------------------

% read gene ids and start and end locations
genelocs=read_stringfile('~/synergy_data/PolII/Mapping_results/genes.txt',[32 9 10 13]);
genelocs={genelocs{2:end}};
n_genes=length(genelocs);
geneids=zeros(n_genes,1);
for k=1:n_genes,
  geneids(k)=str2double(genelocs{k}{1}(6:end-1));
end;

% read gene ids and strand; match them to the previous gene list
genestrands_temp=read_stringfile('~/synergy_data/PolII/Mapping_results/genestrands_no_lrg_genes.txt',[32 9 10 13]);
geneids_temp=zeros(length(genestrands_temp),1);
for k=1:length(genestrands_temp),
  geneids_temp(k)=str2double(genestrands_temp{k}{3}(6:end-1));
end;
genestrands=zeros(n_genes,1);
for k=1:n_genes,
  l=find(geneids_temp==geneids(k));
  if length(l)==0,
    genestrands(k)=NaN;
  else    
    l=l(1);
    genestrands(k)=str2double(genestrands_temp{l}{2});
  end;
end;


%-----------------------------
% Create bin array, one bin for each gene
%-----------------------------
bininfo=zeros(length(genelocs),6);
for i=1:length(genelocs),
  genetemp=genelocs{i};

  % get chromosome index in chromosome list
  genechromosome=genetemp{2}(2:end-1);
  line_chrindex=-1;
  for k=1:n_chromosomes,
    if strcmp(genechromosome,chromosomenames{k})==1,
      line_chrindex=k;
      break;
    end;
  end;

  bininfo(i,1)=line_chrindex; % chromosome index
  bininfo(i,2)=str2double(genetemp{3}(2:end-1)); % start location
  bininfo(i,3)=str2double(genetemp{4}(2:end-1)); % end location
  bininfo(i,4)=i; % line number in an original file, never mind this
  bininfo(i,5)=geneids(i);  % ENSEMBL id of the gene
  bininfo(i,6)=genestrands(i); % strand identifier (+1 or -1)
end;



% Omit genes for which we can't find the strand or the chromosome...
% I=find((isnan(bininfo(:,6))==0) | (bininfo(:,1)==-1));
% bininfo=bininfo(I,:);



%-----------------------------
% Sort bin array in ascending order of start location. 
% The POL2 activity computation code assumes bins have been sorted like that.
%-----------------------------
[y,I1]=sort(bininfo(:,2));
bininfo=bininfo(I1,:);



%-----------------------------
% For each time point, count amount of POL2-read-basepairs that 
% overlap each bin, weighted by scores of the reads.
%-----------------------------


%allbins=zeros(n_genes,10);
filenames={'pol0.mat','pol5.mat','pol10.mat','pol20.mat','pol40.mat','pol80.mat','pol160.mat','pol320.mat','pol640.mat','pol1280.mat'};
dvalues=[190 190 196 192 185 189 201 205 194 189];


% load up all data
allgenebins=cell(size(bininfo,1),10);
for timepoint=1:10,
  load(filenames{timepoint});
  d=dvalues(timepoint);
  max_duplicates=2;  
  subbin_length=200;
  %binsthistime=zeros(n_genes,1);
  
  for chr_index=1:n_chromosomes,
    timepoint
    chr_index
    I=find(bininfo(:,1)==chr_index);
    tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int8(bininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
    for k=1:length(I),
      allgenebins{I(k),timepoint}=tempbins{k};
    end;
  end;
  clear temppol;
end;

save all_gene_pol2bins.mat allgenebins bininfo
