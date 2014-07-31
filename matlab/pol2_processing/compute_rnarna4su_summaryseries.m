beep_on_error(1);
warning("error",'Octave:divide-by-zero');


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
load rnaurna.mat

% perform normalization with eight known stable-to-treatment genes
normalizergenes={'"Ppp2r1a"','"Ndufs5"','"Psma7"','"Tomm7"','"Psmb4"','"Ndufa7"','"Eif4h"','"Capza1"'};
% find normalizer indices
normalizerindices=nan*ones(length(normalizergenes),1);
normalizerseries_rna4su=nan*ones(length(normalizergenes),length(measurementtimes));
normalizerseries_rnatotal=nan*ones(length(normalizergenes),length(measurementtimes));
for i=1:length(normalizergenes),
  for l=1:length(genenames),
    if strcmp(normalizergenes{i},genenames{l})==1,
      normalizerindices(i)=l;
      normalizerseries_rna4su(i,:)=pol_summaryseries(l,:);
      normalizerseries_rnatotal(i,:)=rna_summaryseries(l,:);
    end;
  end;
end;

pol_summaryseries=pol_summaryseries./repmat(mean(normalizerseries_rna4su,1),[size(pol_summaryseries,1) 1]);
rna_summaryseries=rna_summaryseries./repmat(mean(normalizerseries_rnatotal,1),[size(rna_summaryseries,1) 1]);

save rnaurna_normalized.mat pol_summaryseries rna_summaryseries genenames genenanoids measurementtimes -mat



beep();
