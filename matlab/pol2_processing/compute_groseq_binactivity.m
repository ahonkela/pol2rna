mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/'
groseqdir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/GROseq/Mapping_results/'

%mybasedir='~/synergy_data/tempcodebranch/';
%h3k4me3dir='~/synergy_data/H3K4me3/Mapping_results/'
%pol2dir='~/synergy_data/PolII/Mapping_results/'
%groseqdir='~/synergy_data/GROseq/Mapping_results/'


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


if 0,
%-----------------------------
% Code for converting .bed files of GROseq data into matlab cell
% arrays. Only needs to be run once.
%-----------------------------

cd(mexcodedir)
mex read_mappingfile.c

cd(groseqdir)  
gunzipcmd='/opt/local/bin/gunzip '
%gunzipcmd='/bin/gunzip '

system([gunzipcmd '-c GROseq_Vehicle_rep1.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_0a.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

system([gunzipcmd '-c GROseq_Vehicle_rep2.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_0b.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

system([gunzipcmd '-c GROseq_E2_10m_rep1.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_10a.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

system([gunzipcmd '-c GROseq_E2_10m_rep2.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_10b.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

system([gunzipcmd '-c GROseq_E2_40m_rep1.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_40a.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

system([gunzipcmd '-c GROseq_E2_40m_rep2.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_40b.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

system([gunzipcmd '-c GROseq_E2_160m_rep1.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_160a.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

system([gunzipcmd '-c GROseq_E2_160m_rep2.bed.gz > tempfile.bed'])
tempgroseq=read_mappingfile('tempfile.bed',int32(0));
save groseq_160b.mat tempgroseq
system('/bin/rm tempfile.bed')
tempgroseq=1;

end;



if 1,
  cd(mexcodedir)
  mex compute_pol2activityovergenes_c.c

  cd(pol2dir)
  load pol2_for_matti_ver3.mat % provides variable bininfo
  
  
%  chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM','rRNA'};
  chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
  n_chromosomes=length(chromosomenames);
  
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
  allgenebins=cell(size(bininfo,1),ntimepoints);
  fprintf(1,'Initializing allgenebins\n');
  %load all_gene_h3k4me3bins.mat
  %allgenebins=h3k4me3bins;
  
  for timepoint=1:ntimepoints,
  %for timepoint=9:ntimepoints,
    fprintf(1,'Loading data of time point %d\n',timepoint);
    load(filenames{timepoint}); % provides variable tempgroseq

    % sort reads by start position
    for k=1:size(tempgroseq,1),
      [y,I]=sort(tempgroseq{k,1});
      tempgroseq{k,1}=tempgroseq{k,1}(I); % start position in chromosome k
      tempgroseq{k,2}=tempgroseq{k,2}(I); % end position in chromosome
      tempgroseq{k,3}=tempgroseq{k,3}(I); % strand
      tempgroseq{k,4}=tempgroseq{k,4}(I); % score
      tempgroseq{k,7}=tempgroseq{k,7}(I); % original line number in .BED file
    end;
    
    % set read scores to uniform 1 for groseq data
    for k=1:size(tempgroseq,1),
      tempgroseq{k,4}(:)=1;
    end;    
    
    d=dvalues(timepoint);
    max_duplicates=20000000;  
    subbin_length=200;
    %binsthistime=zeros(n_genes,1);
    
    for chr_index=1:n_chromosomes,
      timepoint
      chr_index
      
      % ensure that read scores are in double format, not single      
      if chr_index<=size(tempgroseq,1),
	tempgroseq{chr_index,4}=double(tempgroseq{chr_index,4});
      end;
      
      I=find(bininfo(:,1)==chr_index);
      tempbins=compute_pol2activityovergenes_c(tempgroseq,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int8(bininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allgenebins{I(k),timepoint}=tempbins{k};
      end;
    end;
    clear tempgroseq;
  end;
  groseqbins = allgenebins;
  save all_gene_groseqbins.mat groseqbins -mat
  
end;
