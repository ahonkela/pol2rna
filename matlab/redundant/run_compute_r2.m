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

fid = fopen('~/mlprojects/pol2rnaseq/matlab/file_predictions_config.txt');
conf = textscan(fid, '%s');
fclose(fid);
conf = conf{1};

g = importdata(['~/mlprojects/pol2rnaseq/matlab/', conf{1}]);
myI = mybase:mymod:length(g);

for k=myI,
  fprintf('Running gene %d/%d\n', find(k==myI), length(myI));
  compute_r2_file(g{k}, conf{3}, conf{2});
end
