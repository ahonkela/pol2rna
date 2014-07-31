mybasedir_code='~/mlprojects/';
mybasedir_data='~/projects/pol2rnaseq/';


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


resultdir = [mybasedir_data, 'analyses/hmc_results'];
%d = dir([resultdir, '/rna/*.mat']);
changing_genes = importdata([mybasedir_data, 'data/changing_genes.txt']);
mygenes = changing_genes;

if ~exist('id', 'var'),
  id = '2012-01-23';
end

dic_joint = zeros(size(mygenes));
dic_rna = zeros(size(mygenes));
dic_pol = zeros(size(mygenes));

result_fname = sprintf('~/projects/pol2rnaseq/analyses/dics_%s_%d.mat', id, myI(1));

if exist(result_fname, 'file'),
  load(result_fname);
end

%for k=1:length(mygenes),
for k=myI,
  fprintf('Running gene %d/%d\n', k, length(mygenes));
  if ~any([dic_joint(k), dic_pol(k), dic_rna(k)] == 0),
    continue;
  end
  gene_name = mygenes{k};
  fname = sprintf('%s_samples_%s.mat', gene_name, id);
  try,
    r = load([resultdir, '/joint/', fname]);
    fprintf('loaded %s\n', [resultdir, '/joint/', fname]);
    dic_joint(k) = compute_dic_for_model(r);
    r = load([resultdir, '/pol2/', fname]);
    fprintf('loaded %s\n', [resultdir, '/pol2/', fname]);
    dic_pol(k) = compute_dic_for_model(r);
    r = load([resultdir, '/rna/', fname]);
    fprintf('loaded %s\n', [resultdir, '/rna/', fname]);
    dic_rna(k) = compute_dic_for_model(r);
  catch,
    fprintf('Some results not ready yet.\n')
  end
  fprintf('Joint: %f, Separate: %f\n', dic_joint(k), dic_pol(k) + dic_rna(k));
end

save(result_fname, 'mygenes', 'dic_joint', 'dic_pol', 'dic_rna');
