mybasedir='/Users/hasanogul/jaakkos_files/synergy/mlprojects/';
pol2dir='/Users/hasanogul/jaakkos_files/synergy/synergy_data/PolII/Mapping_results/';

% for kernel-level computations
path1=[mybasedir 'kern/matlab/jaakko_testversion/'];
% for model-level computations
path2=[mybasedir 'gpsim/matlab/jaakko_testversion/'];
% for optimiDefaultConstraint.m
path3=[mybasedir 'optimi/matlab/'];
% for lnDiffErfs.m
path4=[mybasedir 'ndlutil/matlab/'];
% for addPrior.m
path5=[mybasedir 'prior/matlab/'];
% for dist2.m
path6=[mybasedir 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m
path7=[mybasedir 'mltools/matlab/'];
% for various experiment things
path8=[mybasedir 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)

cd(pol2dir)
load pol2_for_matti_ver3.mat
pol_summaryseries=normalizedbygeommeanmedian_pol2series_last20percent;

gene_index=find(bininfo(:,5)==196208);  % GREB1, ENSEMBL-id 196208;

numGenes=1;

times=cell(1,2);

% POL2 observation times
times{1}=[0 5 10 20 40 80 160 320 640 1280]';
%times{1}=[0 1 2 3 4 5 6 7 8 9]';

% mRNA observation times
times{2}=[0 5 10 20 40 80 160 320 640 1280]';
%times{2}=[0 1 2 3 4 5 6 7 8 9]';

dataVals=cell(1,2);
dataVals{1}=pol_summaryseries(gene_index,:)';
dataVals{2}=rna_summaryseries(gene_index,:)';
