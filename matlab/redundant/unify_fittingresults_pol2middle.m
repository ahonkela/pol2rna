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


cd(analysisdir)

allresults_jointmodels={};
allresults_jointtransforminfos={};
allresults_loglikelihoods=[];
allresults_geneindices=[];
allresults_ensemblids=[];
for k=0:10:90,  
  fitname=sprintf('fittingresults_pol2middle_temp_%d.mat',k);
  load(fitname);
  if k==0,
    allresults_jointmodels=results_jointmodels;
    allresults_jointtransforminfos=results_jointtransforminfos;
    allresults_loglikelihoods=results_loglikelihoods;
    allresults_geneindices=results_geneindices;
    allresults_ensemblids=results_ensemblids;
  else
    for l=1:length(results_jointmodels),
      [k l]
      if ~isempty(results_jointmodels{l}),
	if (l>length(allresults_jointmodels)) ...
	      || (isempty(allresults_jointmodels{l})) ...
	      || (allresults_loglikelihoods(l,3)<results_loglikelihoods(l,3)),
	  allresults_jointmodels{l}=results_jointmodels{l};
	  allresults_jointtransforminfos{l}=results_jointtransforminfos{l};
	  allresults_loglikelihoods(l,:)=results_loglikelihoods(l,:);
	  allresults_geneindices(l)=results_geneindices(l);
	  allresults_ensemblids(l)=results_ensemblids(l);
	end;	  
      end;      
    end;    
  end;  
end;



%---------------------------------
% Recompute naive likelihoods
%---------------------------------
for k=1:length(allresults_jointmodels),
  if mod(k,100)==0,
    k
  end;  
    
  if ~isempty(allresults_jointmodels{k}),      
    % Naive log-likelihood of POL2
    datatemp=allresults_jointmodels{k}.y(1:length(allresults_jointmodels{k}.t));
    datatemp=datatemp-min(datatemp);  
    if mean(datatemp.^2)>0,
      datatemp=datatemp/sqrt(mean(datatemp.^2));
      dim = size(datatemp, 1);
      tempK=eye(dim)*var(datatemp);
      tempm=datatemp-mean(datatemp);
      ll1 = -(dim/2)*log(2*pi) - (1/2)*log(det(tempK)) - (1/2)*tempm'*inv(tempK)*tempm;
    else
      ll1=inf;
    end;
    
    % Naive log-likelihood of RNA
    datatemp=allresults_jointmodels{k}.y(length(allresults_jointmodels{k}.t)+1:length(allresults_jointmodels{k}.t)*2);
    if mean(datatemp.^2)>0,
      datatemp=datatemp/sqrt(mean(datatemp.^2));
      dim = size(datatemp, 1);
      tempK=eye(dim)*var(datatemp);
      tempm=datatemp-mean(datatemp);
      ll2 = -(dim/2)*log(2*pi) - (1/2)*log(det(tempK)) - (1/2)*tempm'*inv(tempK)*tempm;
    else
      ll2=inf;
    end;
    
    allresults_loglikelihoods(k,1)=ll1+ll2;
  end;
end;
