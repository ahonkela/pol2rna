%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
%mybasedir_code='~/synergy_data/tempcodebranch/';
%mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
mybasedir_code='~/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/share/work/jtpelto/synergy-data/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';
mybasedir_data='~/projects/pol2rnaseq/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
datadir=[mybasedir_data 'data/'];
analysisdir=[mybasedir_analyses 'analyses/'];

% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/'];
% for optimiDefaultConstraint.m
path3=[mybasedir_code 'optimi/matlab/'];
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



cd(datadir)
%load h3k4me3_series.mat
%load series_for_matti_ver3.mat
load('pol2_for_matti_ver3.mat', 'bininfo');
load('pol2_summaryseries_2012_03.mat');
r = load('rna_new_data.mat');

% gene_index=find(bininfo(:,5)==196208);  % GREB1, ENSEMBL-id 196208;

cd(analysisdir)

[I, A, B] = intersect(r.geneids, bininfo(:, 5));

interestinggenes = A; %find(r.pvals(A) < 0.1);


results_geneindices=[];
results_ensemblids=[];
results_loglikelihoods=[];
results_jointmodels={};
results_jointtransforminfos={};

n_interesting_genes=length(interestinggenes);
n_interesting_genes
%startpercent
%kstart=floor(startpercent*n_interesting_genes/100)+1;
%kend=n_interesting_genes;
%kstart=2;kend=2;
kstart = myindex;
kstep = mycount;
kend = n_interesting_genes;
for k=kstart:kstep:kend,
  rna_index = interestinggenes(k); 
  gene_index= B(interestinggenes(k));
  
  randn('seed',bininfo(gene_index,5));
  rand('seed',bininfo(gene_index,5)+1234567);
  
  dataVals1=pol2_summaryseries(gene_index,:)';
  dataVals2=double(r.counts(rna_index,:)');
  timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
    
  maxiters=100; ninits=10; lengthscale=2;
  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels_celltimes_newdata({temptimes,temptimes},{tempvals1,tempvals2},lengthscale,maxiters,ninits,[],[]);
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

  
  
  if mod((k-kstart)/kstep,5)==0,
    tempfilename=sprintf('fittingresults_shifted_longerdelay_new_%s_%d_%d.mat', id, myindex, mycount);
    safeSave(tempfilename,'results_geneindices','results_ensemblids',...
       'results_loglikelihoods','results_jointmodels','results_jointtransforminfos');
  end;
end;

tempfilename=sprintf('fittingresults_shifted_longerdelay_new_%s_%d_%d.mat', id, myindex, mycount);
safeSave(tempfilename,'results_geneindices','results_ensemblids',...
     'results_loglikelihoods','results_jointmodels','results_jointtransforminfos');
