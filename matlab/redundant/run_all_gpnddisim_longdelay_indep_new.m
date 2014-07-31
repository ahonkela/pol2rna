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

  %------------------------  
  % version with POL2 only
  %------------------------  
  maxiters=100; ninits=8; lengthscale=2;
  %lengthscale=20;
  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels_celltimes_newdata({temptimes,[]},{tempvals1,[]},lengthscale,maxiters,ninits,[],[]);

  fprintf(1, 'Gene %d (ENSG %d) POL2 only: optimized loglik: %f\n',...
	  k,bininfo(gene_index,5),joint_ll);

  results_geneindices(k)=gene_index;
  results_ensemblids(k)=bininfo(gene_index,5);
  results_pol2_loglikelihoods(k,:)=[naive_ll rbf_ll joint_ll];
  results_pol2_jointmodels{k}=jointmodelb;
  results_pol2_jointtransforminfos{k}=jointtransforminfo;



  %------------------------  
  % version with RNA only
  %------------------------  
  maxiters=100; ninits=10; lengthscale=2;
  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels_celltimes_newdata({[],temptimes},{[],tempvals2},lengthscale,maxiters,ninits,[],[]);

  fprintf(1, 'Gene %d (ENSG %d) RNA only: optimized loglik: %f\n',...
	  k,bininfo(gene_index,5),joint_ll);

  results_geneindices(k)=gene_index;
  results_ensemblids(k)=bininfo(gene_index,5);
  results_rna_loglikelihoods(k,:)=[naive_ll rbf_ll joint_ll];
  results_rna_jointmodels{k}=jointmodelb;
  results_rna_jointtransforminfos{k}=jointtransforminfo;

  
  
  if mod((k-kstart)/kstep,5)==0,
    tempfilename=sprintf('fittingresults_shifted_polorrna_new_%s_%d_%d.mat', id, myindex, mycount);
    safeSave(tempfilename,'results_geneindices','results_ensemblids',...
       'results_pol2_loglikelihoods','results_pol2_jointmodels','results_pol2_jointtransforminfos',...
       'results_rna_loglikelihoods','results_rna_jointmodels','results_rna_jointtransforminfos');
  end;
end;

tempfilename=sprintf('fittingresults_shifted_polorrna_new_%s_%d_%d.mat', id, myindex, mycount);
safeSave(tempfilename,'results_geneindices','results_ensemblids',...
   'results_pol2_loglikelihoods','results_pol2_jointmodels','results_pol2_jointtransforminfos',...
   'results_rna_loglikelihoods','results_rna_jointmodels','results_rna_jointtransforminfos');



if 0,
  k=2;
  model = results_jointmodels{k};
  ensemblid=results_ensemblids(k);
  loglik=results_loglikelihoods(k);
  plottitle=makeplottitle(model,loglik,ensemblid);
  predicttimes=timeshift + 1280*(([0:256]'/256).^2);
  plotpredictions(model,predicttimes,2,1,0,1,1,0,[],min(model.t));
  h=gca;
  ylim=get(h,'ylim');
  xlim=get(h,'xlim');
  xloc=xlim(1);
  yloc=ylim(1)+(ylim(2)-ylim(1))*2.6238;
  h=text(xloc,yloc,plottitle);
  set(h,'fontname','Helvetica');
  drawnow;
  print(sprintf('ENSG%011d_GP_shifted_delay%f.eps',ensemblid,model.delay),'-depsc','-S1024,800');
end;

