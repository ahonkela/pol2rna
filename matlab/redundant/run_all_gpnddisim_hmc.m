%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
%mybasedir_code='~/synergy_data/tempcodebranch/';
mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
mybasedir_data='/share/work/jtpelto/synergy-data/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
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




%---------------------------------------------------
% Load maximum likelihood fitting results with small and long delays
%---------------------------------------------------
cd(analysisdir)

load allresults_shifted_polorrna4.mat
% allresults_ensemblids_pol2 allresults_geneindices_pol2 allresults_jointmodels_pol2 allresults_jointtransforminfos_pol2 allresults_loglikelihoods_pol2
% allresults_ensemblids_rna allresults_geneindices_rna allresults_jointmodels_rna allresults_jointtransforminfos_rna allresults_loglikelihoods_rna
allresults_jointtransforminfos_pol2 = allresults_tinfos_pol2;
allresults_jointtransforminfos_rna = allresults_tinfos_rna;


load allresults_shifted_longerdelay5.mat
allresults_ensemblids_joint = allresults_ensemblids; 
allresults_geneindices_joint = allresults_geneindices;
allresults_jointmodels_joint = allresults_jointmodels;
allresults_jointtransforminfos_joint = allresults_jointtransforminfos;
allresults_loglikelihoods_joint = allresults_loglikelihoods;




%---------------------------------------------------
% Compute estimated probability of joint fit being useful
%---------------------------------------------------

probcomparison = nan*ones(max([length(allresults_jointtransforminfos_pol2) ...
     length(allresults_jointtransforminfos_joint)]),1);
for i=1:length(probcomparison),
  if (i<=size(allresults_loglikelihoods_joint,1)) ...
    && (i<=size(allresults_loglikelihoods_pol2,1)),
    if (~isempty(allresults_jointmodels_joint{i})) && (~isempty(allresults_jointmodels_pol2{i})),
      probcomparison(i)=allresults_loglikelihoods_joint(i,3) ...
          -allresults_loglikelihoods_pol2(i,3) -allresults_loglikelihoods_rna(i,3);
    end;
  end;  
end;


% Rank genes by probability of joint fit being useful
[y,I]=sort(-probcomparison);







cd(h3k4me3dir)
%load h3k4me3_series.mat
%load series_for_matti_ver3.mat
load series_for_matti_corrected_ver1.mat


cd(analysisdir)
load models_Pol2.mat
templls=lls_Pol2(:,1)-lls_Pol2(:,3);
[y,I]=sort(-templls);
interestinggenes=genes_Pol2;
interestinggenes=interestinggenes(I);
mattisparameters=params_Pol2;
mattisparameters=mattisparameters(I,:);

mattisparameters_all=nan*ones(size(bininfo,1),size(mattisparameters,2));
mattisparameters_all(interestinggenes,:)=mattisparameters;


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
end;  


if 0,
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
  %interesting_ensemblids=[229732 142207 140961 185818 196208 211891 107672];
  %interesting_ensemblids=[196208 211891 107562];
  %interesting_ensemblids=[229732 142207 140961 185818];
  interesting_ensemblids=[196208];
  %interesting_ensemblids=[229732 142207];
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
n_interesting_genes
%startpercent
kstart=floor(startpercent*n_interesting_genes/100)+1;
kend=n_interesting_genes;
%kstart=2;kend=2;
for k=kstart:kend,
  gene_index=interestinggenes(k);
  
  randn('seed',bininfo(gene_index,5));
  rand('seed',bininfo(gene_index,5)+1234567);
  
  dataVals1=pol_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,:)';
  timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
    
  maxiters=100; ninits=8; lengthscale=2;

  initializationtype=1;
  parameterranges=[];
  use_fixedrnavariance=[];
  [jointmodel,temptransforminfo]=createNdSimDisim_longdelay(timevector,dataVals1,dataVals2,lengthscale,initializationtype,parameterranges,use_fixedrnavariance);

  samples = gpnddisimSampleHMC(jointmodel, 2, 1000);

  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels_longdelay(temptimes,tempvals1,tempvals2,lengthscale,maxiters,ninits,[],[]);
  %  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createGeneGPModelsConditional(timevector,dataVals1,dataVals2,lengthscale,maxiters,ninits);

  % try a run from Matti's parameters
  [pars1,nams1]=gpnddisimExtractParam(jointmodelb)
  parst1=transformParametersWithSettings(pars1,jointtransforminfo,'atox');  
  parst2=parst1;
  parst2(9)=mattisparameters_all(gene_index,1); % RNA start mean
  parst2(8)=mattisparameters_all(gene_index,2); % Basal rate
  parst2(4)=mattisparameters_all(gene_index,3)^2; % DISIM-level variance = square of sensitivity
  parst2(3)=mattisparameters_all(gene_index,4); % DISIM-level decay
  parst2(5)=mattisparameters_all(gene_index,5); % DISIM-level delay
  parst2(10)=jointmodelb.y(1); % POL2 start mean
  parst2(6)=1e-3; % start with low POL2 noise assumption
  parst2(7)=1e-3; % start with low RNA noise assumption
  pars2=transformParametersWithSettings(parst2,jointtransforminfo,'xtoa');
  model2=gpnddisimExpandParam(jointmodelb,pars2);
  model2_ll=gpnddisimLogLikelihood(model2);
  
  % plotpredictions(model2,[0:5:1280]',2,1,1,'exampletitle');
  % drawnow;
  
  try  
    model3=gpnddisimOptimise(model2,1,maxiters);
  catch
    model3=model2;
  end;
  model3_ll=gpnddisimLogLikelihood(model3);
  if (model3_ll > joint_ll),
    jointmodelb=model3;
    joint_ll=model3_ll;
  end;  
  % plotpredictions(jointmodelb,[0:5:1280]',2,1,1,'exampletitle');
  % drawnow;
  fprintf(1, 'Gene %d (ENSG %d): optimized loglik: %f, ll w. Mattis params %f, ll opt w. M.params %f\n',...
	  k,bininfo(gene_index,5),joint_ll,model2_ll,model3_ll);
  %pause

  % plotpredictions(results_jointmodels{k},[0:5:1280]',2,1,1,sprintf('ENSG %d',results_ensemblids(k)));

  results_geneindices(k)=gene_index;
  results_ensemblids(k)=bininfo(gene_index,5);
  results_loglikelihoods(k,:)=[naive_ll rbf_ll joint_ll];
  results_jointmodels{k}=jointmodelb;
  results_jointtransforminfos{k}=jointtransforminfo;

  
  
  if mod(k,5)==0,
    tempfilename=sprintf('fittingresults_shifted_longerdelay_startpercent%d.mat',startpercent);
    save(tempfilename,'results_geneindices','results_ensemblids',...
       'results_loglikelihoods','results_jointmodels','results_jointtransforminfos','-mat');
  end;
end;

tempfilename=sprintf('fittingresults_shifted_longerdelay_startpercent%d.mat',startpercent);
save(tempfilename,'results_geneindices','results_ensemblids',...
     'results_loglikelihoods','results_jointmodels','results_jointtransforminfos','-mat');



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

