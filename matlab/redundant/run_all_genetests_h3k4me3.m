%mybasedir_code='/users/hasanogul/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
mybasedir_code='~/synergy_data/tempcodebranch/';

%mybasedir_data='/users/hasanogul/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
mybasedir_data='~/synergy_data/';

% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/jaakko_testversion/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/jaakko_testversion/'];
% for optimiDefaultConstraint.m
path3=[mybasedir_code 'optimi/matlab/jaakko_testversion/'];
% for lnDiffErfs.m
path4=[mybasedir_code 'ndlutil/matlab/'];
% for addPrior.m
path5=[mybasedir_code 'prior/matlab/'];
% for dist2.m
path6=[mybasedir_code 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m
path7=[mybasedir_code 'mltools/matlab/'];
% for various experiment things
path8=[mybasedir_code 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)


%cd(pol2dir)
%load pol2_for_matti_ver3.mat
%pol_summaryseries=normalizedbygeommeanmedian_pol2series_last20percent;
cd(h3k4me3dir)
%load h3k4me3_series.mat
load series_for_matti.mat

% gene_index=find(bininfo(:,5)==196208);  % GREB1, ENSEMBL-id 196208;


if 0,
%---------------------------------------------------
% Find genes with: 
% - substantial POL2 activity 
% - substantial POL2 variance over time
% - substantial RNA activity 
% - substantial RNA variance over time
%---------------------------------------------------

pol_summaryseries_means=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
  pol_summaryseries_means(k)=mean(pol_summaryseries(k,:));
end;
pol_summaryseries_variances=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
  pol_summaryseries_variances(k)=var(pol_summaryseries(k,:));
end;
rna_summaryseries_means=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
  rna_summaryseries_means(k)=mean(rna_summaryseries(k,:));
end;
rna_summaryseries_variances=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
  rna_summaryseries_variances(k)=var(rna_summaryseries(k,:));
end;

Interestinggenes=find((pol_summaryseries_means>=exp(9)) & ...
                      (pol_summaryseries_variances>=exp(15)) & ...
                      (rna_filteringresults==1));


%------------------
% Alternative selection: use genes with early RNA peaks as interesting
%------------------
analysisdir='/Users/hasanogul/jaakkos_files/synergy/analyses/';
stringfile=read_stringfile([analysisdir 'earlygenes_t5.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t5=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t10.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t10=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t20.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t20=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t40.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t40=tempgenes;
stringfile=read_stringfile([analysisdir 'earlygenes_t80.txt'], [32 9 10 13]);
tempgenes=zeros(length(stringfile)-1,1);
for k=2:length(stringfile), tempgenes(k-1)=str2double(stringfile{k}{1}(6:(end-1))); end;
interestinggenes_t80=tempgenes;
interesting_ensemblids=unique([interestinggenes_t5;interestinggenes_t10;interestinggenes_t20;interestinggenes_t40;interestinggenes_t80]);

end;


if 0,
interesting_ensemblids=[104267 105219 112658 196208];

interestinggenes=[];
for k=1:length(interesting_ensemblids),
  gene_index=find(bininfo(:,5)==interesting_ensemblids(k));
  if ~isempty(gene_index),
    interestinggenes=[interestinggenes gene_index(1)];
  end;
end;
end;






results_geneindices=[];
results_ensemblids=[];
results_loglikelihoods=[];
results_jointmodels={};
results_jointtransforminfos={};

n_interesting_genes=length(interestinggenes);
kstart=1;
kend=n_interesting_genes;
for k=kstart:kend,
  gene_index=interestinggenes(k);
  
  randn('seed',bininfo(gene_index,5));
  
  %dataVals1=pol_summaryseries(gene_index,:)';
  %dataVals2=rna_summaryseries(gene_index,:)';
  %timevector=[0 5 10 20 40 80 160 320 640 1280]';
  dataVals1=h3k4me3_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,[1 3 4 5 6 7 8 9 10])';
  timevector=[0 10 20 40 80 160 320 640 1280]';

  dataVals1=dataVals1(1:7);
  dataVals2=dataVals2(1:7);
  timevector=timevector(1:7);
  
  maxiters=100; ninits=1; lengthscale=5;
  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createGeneGPModels(timevector,dataVals1,dataVals2,lengthscale,maxiters,ninits);
%  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createGeneGPModelsConditional(timevector,dataVals1,dataVals2,lengthscale,maxiters,ninits);
  
  results_geneindices(k)=gene_index;
  results_ensemblids(k)=bininfo(gene_index,5);
  results_loglikelihoods(k,:)=[naive_ll rbf_ll joint_ll];
  results_jointmodels{k}=jointmodelb;
  results_jointtransforminfos{k}=jointtransforminfo;
end;


