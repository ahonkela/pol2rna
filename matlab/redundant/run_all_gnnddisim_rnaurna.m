%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
mybasedir_code='~/synergy_data/tempcodebranch/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
mybasedir_data='~/synergy_data/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
analysisdir=[mybasedir_analyses 'analyses/'];
rnaurnadir=[mybasedir_data 'RNAURNA/'];

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



cd(rnaurnadir)
%load h3k4me3_series.mat
%load series_for_matti_ver3.mat
load rnaurna.mat

interestinggenes=[1:size(pol_summaryseries,1)];






results_geneindices=[];
results_ensemblids=[];
results_loglikelihoods=[];
results_jointmodels={};
results_jointtransforminfos={};

n_interesting_genes=length(interestinggenes);
n_interesting_genes
startpercent
kstart=floor(startpercent*n_interesting_genes/100)+1;
kend=n_interesting_genes;
%kstart=2;kend=2;
for k=kstart:kend,
  gene_index=interestinggenes(k);
  
  randn('seed',bininfo(gene_index,5));
  rand('seed',bininfo(gene_index,5)+1234567);
  
  dataVals1=pol_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,:)';
  timevector=measurementtimes; %[0 5 10 20 40 80 160 320 640 1280]';
  %dataVals1=h3k4me3_summaryseries(gene_index,:)';
  %dataVals2=rna_summaryseries(gene_index,[1 3 4 5 6 7 8 9 10])';
  %timevector=[0 10 20 40 80 160 320 640 1280]';

  % dataVals1=dataVals1(1:8);
  % dataVals2=dataVals2(1:8);
  % timevector=timevector(1:8);

  % temptimes=timevector([1 (3:length(timevector))]);
  % tempvals1=dataVals1([1 (3:length(timevector))]);
  % tempvals2=dataVals2([1 (3:length(timevector))]);  

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;

  % normalization
  % tempvals1=tempvals1-min(tempvals1);
  % if mean(tempvals1.^2)>0,
  %   tempvals1=tempvals1/sqrt(mean(tempvals1.^2));
  % end;
  % if mean(tempvals2.^2)>0,
  %   tempvals2=tempvals2/sqrt(mean(tempvals2.^2));
  % end;
    
  maxiters=100; ninits=8; lengthscale=5;
  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels(temptimes,tempvals1,tempvals2,lengthscale,maxiters,ninits);
  %  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createGeneGPModelsConditional(timevector,dataVals1,dataVals2,lengthscale,maxiters,ninits);

  % plotpredictions(jointmodelb,[0:5:1280]',2,1,1,'exampletitle');
  % drawnow;
  fprintf(1, 'Gene %d (ENSG %d): optimized loglik: %f\n',...
	  k,bininfo(gene_index,5),joint_ll);
  %pause

  % plotpredictions(results_jointmodels{k},[0:5:1280]',2,1,1,sprintf('ENSG %d',results_ensemblids(k)));

  results_geneindices(k)=gene_index;
  results_ensemblids(k)=bininfo(gene_index,5);
  results_loglikelihoods(k,:)=[naive_ll rbf_ll joint_ll];
  results_jointmodels{k}=jointmodelb;
  results_jointtransforminfos{k}=jointtransforminfo;
  if mod(k,5)==0,
    tempfilename=sprintf('fittingresults_rnaurna_temp_%d.mat',startpercent);
    save(tempfilename,'results_geneindices','results_ensemblids',...
       'results_loglikelihoods','results_jointmodels','results_jointtransforminfos','-mat');
  end;
end;



