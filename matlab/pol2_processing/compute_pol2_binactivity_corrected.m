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


if 0,
%-----------------------------
% Code for converting .bed files of POL2 data into matlab cell
% arrays. Only needs to be run once.
%-----------------------------

cd(mexcodedir)
mex read_mappingfile.c

cd(pol2dir)  
gunzipcmd='/opt/local/bin/gunzip '
rmcmd='/bin/rm '

system([gunzipcmd '-c MCF7_No_Treat_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0)); % change 0 to
                                                   % the number of
                                                   % header lines
                                                   % in the .bed file!
save pol0.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_5min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol5.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_10min_PolII_rep_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol10.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_20min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol20.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_40min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol40.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_80min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol80.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_160min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol160.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_320min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol320.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_640min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol640.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

system([gunzipcmd '-c MCF7_E2_1280min_PolII_unique.bed.gz > tempfile.bed'])
temppol=read_mappingfile('tempfile.bed',int32(0));
save pol1280.mat temppol -mat
system([rmcmd 'tempfile.bed'])
temppol=1;

end;



if 1,
  cd(mexcodedir)
  mex compute_pol2activityovergenes_c.c

  cd(pol2dir)
  load pol2_for_matti_ver3.mat   % provides variable bininfo
  
  
  chromosomenames={'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY','chrM'};
  n_chromosomes=length(chromosomenames);
  
  cd(pol2dir)
  
  
  %-----------------------------
  % For each time point, count amount of POL2-read-basepairs that 
  % overlap each bin, weighted by scores of the reads.
  %-----------------------------
  
  %allbins=zeros(n_genes,10);
  filenames={'pol0.mat','pol5.mat','pol10.mat','pol20.mat','pol40.mat','pol80.mat','pol160.mat','pol320.mat','pol640.mat','pol1280.mat'};
  dvalues=[190 190 196 192 185 189 201 205 194 189]; % dvalues obtained from MACS output
  
  ntimepoints=length(filenames);
  % load up all data
  allgenebins=cell(size(bininfo,1),ntimepoints);
  fprintf(1,'Initializing allgenebins\n');
%  load all_gene_h3k4me3bins.mat
%  allgenebins=h3k4me3bins;
  
  for timepoint=1:ntimepoints,
  %for timepoint=9:ntimepoints,
    fprintf(1,'Loading data of time point %d\n',timepoint);
    load(filenames{timepoint});  % provides variable temppol
    d=dvalues(timepoint);
    max_duplicates=2;  
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
      
      I=find(bininfo(:,1)==chr_index);
      tempbins=compute_pol2activityovergenes_c(temppol,int32(chr_index),int32(length(I)),int32(bininfo(I,2)),int32(bininfo(I,3)),int8(bininfo(I,6)),int32(d),int32(max_duplicates),int32(subbin_length));
      for k=1:length(I),
        allgenebins{I(k),timepoint}=tempbins{k};
      end;
    end;
    clear temppol;
  end;
  pol2bins = allgenebins;
  save all_gene_pol2bins.mat pol2bins -mat
  
end;
