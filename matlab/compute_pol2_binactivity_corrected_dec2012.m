mybasedir='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%h3k4me3dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/H3K4me3/Mapping_results/'
pol2gzdir = '/share/synergy/data/2012-03_PolII/Mapping_results/';
pol2dir='/share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed/'



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

mexcodedir=path8


if 0,
%-----------------------------
% Code for converting .bed files of POL2 data into matlab cell
% arrays. Only needs to be run once.
%-----------------------------

cd(mexcodedir)
mex read_mappingfile.c

cd(pol2dir)  
gunzipcmd='/bin/gunzip '
rmcmd='/bin/rm '

if 0,
gzname='MCF7_NoTreat_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol0_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

gzname='MCF7_5min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol5_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

gzname='MCF7_10min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol10_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

gzname='MCF7_20min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol20_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;
end;

gzname='MCF7_40min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol40_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
clear temppol;

gzname='MCF7_80min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol80_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
clear temppol;

gzname='MCF7_160min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol160_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
clear temppol;

gzname='MCF7_320min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol320_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
clear temppol;

gzname='MCF7_640min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol640_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
clear temppol;

gzname='MCF7_1280min_PolII_2012-03_unique.bed.gz';
systemcommand=sprintf('%s -c %s%s > %stempfile.bed',gunzipcmd,pol2gzdir,gzname,pol2dir);
system(systemcommand);
temppol=read_mappingfile('tempfile.bed',int32(4)); % change 0 to the number of header lines in the .bed file!
save pol1280_2012_03.mat temppol -mat
system([rmcmd 'tempfile.bed'])
clear temppol;

end;



if 1,
  cd(mexcodedir)
  mex compute_pol2activityovergenes_c.c

  cd(pol2dir)

  % provides variable bininfo, corrected with active transcript area for some genes
  load /share/mi/workspace/jtpelto/synergy/synergy_data/analyses/new_bininfo/bininfo_dec2012.mat

  % Update part of the bininfo structure based on transcript activities...
  % Compute also a reduced bininfo structure containing only the problem genes
  problems=load('/share/mi/workspace/jtpelto/synergy/synergy_data/analyses/problemgenes.mat');
  % bininfo_problemgenes = nan*ones(problems.nproblemgenes,size(bininfo,2));
  for k=1:problems.nproblemgenes,
    geneid=problems.problemgenes(k);
    bininfoid=find(bininfo(:,5)==geneid);
    bininfo(bininfoid,2)=problems.problemgene_extendedstarts(k);
    bininfo(bininfoid,3)=problems.problemgene_extendedends(k);
    % bininfo_problemgenes(k,:) = bininfo(bininfoid,:);
  end;

  % Sort the reduced bininfo structure by the start position index, required by the C code!
  bininfo_problemgenes = bininfo;
  [ytempsort,Itempsort] = sort(bininfo_problemgenes(:,2));
  bininfo_problemgenes = bininfo_problemgenes(Itempsort,:);
  
  
  chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
  n_chromosomes=length(chromosomenames);
  
  cd(pol2dir)
  
  
  %-----------------------------
  % For each time point, count amount of POL2-read-basepairs that 
  % overlap each bin, weighted by scores of the reads.
  %-----------------------------
  
  %allbins=zeros(n_genes,10);
  filenames={'pol0_2012_03.mat','pol5_2012_03.mat','pol10_2012_03.mat','pol20_2012_03.mat','pol40_2012_03.mat','pol80_2012_03.mat','pol160_2012_03.mat','pol320_2012_03.mat','pol640_2012_03.mat','pol1280_2012_03.mat'};

  % dvalues=[230 214 219 213 217 223 215 208 213 212]; % dvalues obtained from MACS output
  dvalues=[0 0 0 0 0 0 0 0 0 0]; % no shifting or lengthening of reads
  
  ntimepoints=length(filenames);
  % load up all data
  allgenebins_problemgenes=cell(size(bininfo_problemgenes,1),ntimepoints);
  fprintf(1,'Initializing allgenebins_problemgenes\n');
%  load all_gene_h3k4me3bins.mat
%  allgenebins=h3k4me3bins;
  
  for timepoint=1:ntimepoints,
  %for timepoint=9:ntimepoints,
    fprintf(1,'Loading data of time point %d\n',timepoint);
    load(filenames{timepoint});  % provides variable temppol
    d=dvalues(timepoint);

    % For the 2012-03 data, MACS seems to keep only 1 duplicate; 
    % let's do the same even if we do not use MACS for shifting/lengthening reads
    max_duplicates=1;  

    subbin_length=200;
    %binsthistime=zeros(n_genes,1);
    
    for chr_index=1:n_chromosomes,
%    for chr_index=22:n_chromosomes,
      timepoint
      chr_index
      
      % ensure that read scores are in double format, not single
      if chr_index<=size(temppol,1),
        temppol{chr_index,4}=double(temppol{chr_index,4});
      end;
      
      I=find(bininfo_problemgenes(:,1)==chr_index);
      tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo_problemgenes(I,2)),int32(bininfo_problemgenes(I,3)),int8(bininfo_problemgenes(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allgenebins_problemgenes{I(k),timepoint}=tempbins{k};
      end;
    end;
    clear temppol;
  end;
  pol2bins = allgenebins_problemgenes;

  save all_gene_pol2bins_2013_01_02.mat pol2bins -mat
  save /share/mi/workspace/jtpelto/synergy/synergy_data/analyses/new_bininfo/bininfo_dec2012_corrected.mat bininfo -mat 
end;
