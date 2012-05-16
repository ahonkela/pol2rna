N = 200;
filestem = '/share/synergy/analyses/2012-05-07/fittingresults_shifted_longerdelay_new_2012-05-07b_%d_200.mat';

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

addpath(path1,path2,path3,path4,path5,path6,path7,path8);

load('~/projects/pol2rnaseq/data/pol2_for_matti_ver3.mat', 'bininfo');
r = load('~/projects/pol2rnaseq/data/rna_new_data.mat');
[I, A, B] = intersect(r.geneids, bininfo(:, 5));
symbols = r.symbols(A);

figdir = '/share/synergy/analyses/2012-05-07/figures/';

for k=myI,
    mkdir(figdir, sprintf('%d', k));
    cd(sprintf('%s/%d', figdir, k));
    fprintf('%d/%d\n', find(myI==k), length(myI));
    res = load(sprintf(filestem, k));
    I = (res.results_geneindices ~= 0);
    for l=find(I),
        plot_model(res.results_jointmodels{l}, res.results_ensemblids(l), res.results_loglikelihoods(l), symbols{l});
    end
end
