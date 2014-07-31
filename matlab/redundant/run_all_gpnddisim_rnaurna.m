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
load rnaurna_normalized.mat

interestinggenes=[1:size(pol_summaryseries,1)];






results_geneindices=[];
results_genenames={};
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
  randn('seed',gene_index);
  rand('seed',gene_index+1234567);
  
  dataVals1=pol_summaryseries(gene_index,:)';
  dataVals2=rna_summaryseries(gene_index,:)';
  timevector=measurementtimes'; 
  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
  
  maxiters=100; ninits=8; lengthscale=0.25;
    
  % set parameter ranges
  inversewidth_range=1./([1280 0.1].^2);
  pol2effectvar_range=[0.0002 20];%*var(tempvals1);
  if pol2effectvar_range(2)==pol2effectvar_range(1),
    pol2effectvar_range(2)=pol2effectvar_range(1)+1;
  end;  
  pol2noisevar_range=[0.005 5];%*var(tempvals1);
  if pol2noisevar_range(2)==pol2noisevar_range(1),
    pol2noisevar_range(2)=pol2noisevar_range(1)+1;
  end;  
  rnadecay_range=[0.00001 1.5];
  rnaeffectvar_range=[0 20];
  rnadelay_range=[0 5];
  rnanoisevar_range=[0.005 5];%*var(tempvals2);
  if rnanoisevar_range(2)==rnanoisevar_range(1),
    rnanoisevar_range(2)=rnanoisevar_range(1)+1;
  end;
  rnabasal_range=[0 10000];
  rnastartmean_range=[0 10000];
  pol2mean_range=[0 10000];
  parameterranges=[inversewidth_range;pol2effectvar_range;pol2noisevar_range;rnadecay_range;rnaeffectvar_range;rnadelay_range;rnanoisevar_range;rnabasal_range;rnastartmean_range;pol2mean_range];  
      
  use_fixedrnavariance=0;
  [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels(temptimes,tempvals1,tempvals2,lengthscale,maxiters,ninits,parameterranges,use_fixedrnavariance);

  fprintf(1, 'Gene %d (name %s, nanoid %s): optimized loglik: %f\n',...
	  k,genenames{interestinggenes(k)},genenanoids{interestinggenes(k)},joint_ll);

  results_geneindices(k)=gene_index;
  results_genenames{k}=genenames{gene_index};
  results_loglikelihoods(k,:)=[naive_ll rbf_ll joint_ll];
  results_jointmodels{k}=jointmodelb;
  results_jointtransforminfos{k}=jointtransforminfo;
  if mod(k,5)==0,
    tempfilename=sprintf('fittingresults_rnaurna_temp_%d.mat',startpercent);
    save(tempfilename,'results_geneindices','results_genenames',...
       'results_loglikelihoods','results_jointmodels','results_jointtransforminfos','-mat');
  end;
end;

