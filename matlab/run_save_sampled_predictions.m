mybasedir_code='~/mlprojects/';

% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/'];
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

g = importdata('~/mlprojects/pol2rnaseq/matlab/finished_genes_2012-11-13.txt');
myI = mybase:mymod:length(g);

for k=myI,
  fprintf('Running gene %d/%d\n', find(k==myI), length(myI));
  save_sampled_predictions_file(g{k});
end
