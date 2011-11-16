% Analyse how high is a list of genes provided by Stuart Wilson, 
% whose website is at http://www.shef.ac.uk/mbb/staff/wilson
% The list of genes is currently confidential.
% Idea: for each of Stuart's genes, plot a position in a 2D scatterplot,
%   where the coordinates are: the rank of the gene in our "usefulness 
%   of joint GP modeling" ranking, and the Pol2-RNA delay estimated by 
%   the joint GP fit.
% If our estimated delay for Stuart's genes is typically high, and 
% conversely if many of our genes appear in one of Stuart's lists,
% this can indicate that Stuart's explanation for Pol2-RNA delay could
% apply to our data set.
%
% Jaakko Peltonen, Nov 2, 2011.

%mybasedir_code='/media/JPELTONEN4/mlprojects/';
mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';
%mybasedir_code='/share/work/jtpelto/tempsynergy/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';
mybasedir_data='/share/work/jtpelto/synergy-data/';

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
% Load fitting results with small and long delays
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
% Compute estimated probability of gene having long delay
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

[y,Iprobcomparison]=sort(-probcomparison);

