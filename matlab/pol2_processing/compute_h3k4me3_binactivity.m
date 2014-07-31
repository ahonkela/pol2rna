mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/'


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
%-----------------------------
% Code for converting .bed files of POL2 data into matlab cell
% arrays. Only needs to be run once.
%-----------------------------

cd(mexcodedir)
mex read_mappingfile.c

cd(h3k4me3dir)  

system('/opt/local/bin/gunzip -c MCF7_NoTreat_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_0.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_10min_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_10.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_20min_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_20.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_40min_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_40.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_80min_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_80.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_160min_H3K4me_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_160.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_320min_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_320.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_640min_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_640.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

system('/opt/local/bin/gunzip -c MCF7_E2_1280min_H3K4me3_unique.bed.gz > tempfile.bed')
temph3k4me3=read_mappingfile('tempfile.bed',int32(0));
save h3k4me3_1280.mat temph3k4me3 -mat
system('/bin/rm tempfile.bed')
temph3k4me3=1;

end;



if 1,
  cd(mexcodedir)
  mex compute_pol2activityovergenes_c.c

  cd(pol2dir)
  load pol2_for_matti_ver3.mat   % provides variable bininfo
  
  
  chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
  n_chromosomes=length(chromosomenames);
  
  cd(h3k4me3dir)
  
  
  %-----------------------------
  % For each time point, count amount of POL2-read-basepairs that 
  % overlap each bin, weighted by scores of the reads.
  %-----------------------------
  
  
  %allbins=zeros(n_genes,10);
  filenames={'h3k4me3_0.mat','h3k4me3_10.mat','h3k4me3_20.mat','h3k4me3_40.mat','h3k4me3_80.mat','h3k4me3_160.mat','h3k4me3_320.mat','h3k4me3_640.mat','h3k4me3_1280.mat'};
  dvalues=[192 192 196 196 180 197 186 198 211]; % dvalues obtained from MACS output
  
  ntimepoints=length(filenames);
  % load up all data
  allgenebins=cell(size(bininfo,1),ntimepoints);
  fprintf(1,'Initializing allgenebins\n');
%  load all_gene_h3k4me3bins.mat
%  allgenebins=h3k4me3bins;
  
  for timepoint=1:ntimepoints,
  %for timepoint=9:ntimepoints,
    fprintf(1,'Loading data of time point %d\n',timepoint);
    load(filenames{timepoint});  % provides variable temph3k4me3
    d=dvalues(timepoint);
    max_duplicates=2;  
    subbin_length=200;
    %binsthistime=zeros(n_genes,1);
    
    for chr_index=1:n_chromosomes,
      timepoint
      chr_index
      I=find(bininfo(:,1)==chr_index);
      tempbins=compute_pol2activityovergenes_c(temph3k4me3,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int8(bininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allgenebins{I(k),timepoint}=tempbins{k};
      end;
    end;
    clear temph3k4me3;
  end;
  h3k4me3bins = allgenebins;
%  save all_gene_h3k4me3bins.mat h3k4me3bins -mat
  
end;
