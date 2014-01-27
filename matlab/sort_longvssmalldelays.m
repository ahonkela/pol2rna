%mybasedir_code='/media/JPELTONEN4/mlprojects/';
%mybasedir_code='~/synergy_data_desktopbackup/tempcodebranch/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';
mybasedir_code='~/code/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data_desktopbackup/';
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


%---------------------------------------------------
% Load gene data
%---------------------------------------------------
cd(h3k4me3dir)
load series_for_matti_corrected_ver1.mat  % provides pol_summaryseries,rna_summaryseries

cd(analysisdir)
load models_Pol2.mat % provides interestinggenes
templls=lls_Pol2(:,1)-lls_Pol2(:,3);
[y,I]=sort(-templls);
interestinggenes=genes_Pol2;
interestinggenes=interestinggenes(I);


%---------------------------------------------------
% Load fitting results with small and long delays
%---------------------------------------------------
cd(analysisdir)
load allresults_smalldelay.mat

allresults_ensemblids_small = allresults_ensemblids; 
allresults_geneindices_small = allresults_geneindices;
allresults_jointmodels_small = allresults_jointmodels;
allresults_jointtransforminfos_small = allresults_jointtransforminfos;
allresults_loglikelihoods_small = allresults_loglikelihoods;

load allresults_longdelay_aug5.mat

allresults_ensemblids_long = allresults_ensemblids; 
allresults_geneindices_long = allresults_geneindices;
allresults_jointmodels_long = allresults_jointmodels;
allresults_jointtransforminfos_long = allresults_jointtransforminfos;
allresults_loglikelihoods_long = allresults_loglikelihoods;


%---------------------------------------------------
% Compute estimated probability of gene having long delay
%---------------------------------------------------

longdelayprobs = nan*ones(max([length(allresults_jointtransforminfos_small) ...
     length(allresults_jointtransforminfos_long)]),1);
for i=1:length(longdelayprobs),
  if (i<=size(allresults_loglikelihoods_small,1)) ...
    && (i<=size(allresults_loglikelihoods_long,1)),
    temp1=min([allresults_loglikelihoods_long(i,3) ...
               allresults_loglikelihoods_small(i,3)]);
    if (~isempty(allresults_jointmodels_small{i})) && (~isempty(allresults_jointmodels_long{i})),
      longdelayprobs(i)=exp(allresults_loglikelihoods_long(i,3)-temp1) / ...
          (exp(allresults_loglikelihoods_long(i,3)-temp1)+exp(allresults_loglikelihoods_small(i,3)-temp1));    
    end;
  end;  
end;


%---------------------------------------------------
% Compute estimated probability of gene having an effect 
% (compare modeling with long delay to modeling with noise only)
%---------------------------------------------------

effectprobs = nan*ones(max([length(allresults_jointtransforminfos_small) ...
     length(allresults_jointtransforminfos_long)]),1);
for i=1:length(effectprobs),
  if (i<=size(allresults_loglikelihoods_small,1)) ...
    && (i<=size(allresults_loglikelihoods_long,1)),
    temp1=min([allresults_loglikelihoods_long(i,3) ...
               allresults_loglikelihoods_long(i,1)]);
    if (~isempty(allresults_jointmodels_small{i})) && (~isempty(allresults_jointmodels_long{i})),
      effectprobs(i)=exp(allresults_loglikelihoods_long(i,3)-temp1) / ...
          (exp(allresults_loglikelihoods_long(i,3)-temp1)+exp(allresults_loglikelihoods_long(i,1)-temp1));    
    end;
  end;  
end;


%---------------------------------------------------
% Sort score: probability of gene having an effect, 
% times the probability of the effect having long delay.
%---------------------------------------------------
sortscore=longdelayprobs.*effectprobs;
[y,I]=sort(-sortscore);



%---------------------------------------------------
% Put the sorted results in their own structure and save them.
%---------------------------------------------------
sortedresults_loglikelihoods=[];
sortedresults_ensemblids=[];
sortedresults_jointmodels_longdelay={};
sortedresults_jointmodels_smalldelay={};
for k=1:length(I),
  if isnan(sortscore(I(k)))==0,
    m=I(k);
    sortedresults_loglikelihoods(k,1)=allresults_loglikelihoods_long(m,1);
    sortedresults_loglikelihoods(k,2)=allresults_loglikelihoods_small(m,3);
    sortedresults_loglikelihoods(k,3)=allresults_loglikelihoods_long(m,3);
    sortedresults_loglikelihoods(k,4)=sortscore(m);
    sortedresults_loglikelihoods(k,5)=allresults_ensemblids_long(m);
    sortedresults_ensemblids(k)=allresults_ensemblids_long(m);
    sortedresults_jointmodels_longdelay{k}=allresults_jointmodels_long{m};
    sortedresults_jointmodels_smalldelay{k}=allresults_jointmodels_small{m};
    sortedresults_jointtransforminfos_longdelay{k}=allresults_jointtransforminfos_long{m};
    sortedresults_jointtransforminfos_smalldelay{k}=allresults_jointtransforminfos_small{m};
  end;
end;


save smallvslargedelays_ascii.txt sortedresults_loglikelihoods -ascii
save smallvslargedelays.mat sortedresults_loglikelihoods sortedresults_ensemblids ...
   sortedresults_jointmodels_longdelay sortedresults_jointmodels_smalldelay ...
   sortedresults_jointtransforminfos_longdelay sortedresults_jointtransforminfos_smalldelay -mat

