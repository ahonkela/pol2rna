
% N = 200;
% filestem = '~/projects/pol2rnaseq/analyses/fittingresults_shifted_longerdelay_new_2012-05-07b_%d_200.mat';
% filestem = '/share/synergy/analyses/2012-05-07/fittingresults_shifted_longerdelay_new_2012-05-07b_%d_200.mat';

N = 200;
filestem = '/share/synergy/analyses/2012-05-07/fittingresults_shifted_longerdelay_new_2012-05-07b_%d_200.mat';


mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/share/work/jtpelto/synergy-data/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';
% mybasedir_data='~/projects/pol2rnaseq/';

% mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
% datadir=[mybasedir_data 'data/'];
% analysisdir=[mybasedir_analyses 'analyses/'];

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



% if 0,
% load('/share/synergy/data_summaries/pol2_for_matti_ver3.mat', 'bininfo');

figdir = '/share/synergy/analyses/2012-05-07/figures/';



% r = load('~/projects/pol2rnaseq/data/rna_new_data.mat');

rnafile = read_stringfile('/share/synergy/data_summaries/new_rna_counts.txt',[],[' ' 10 13 0]);
%rna_summaryseries=nan*ones(size(bininfo,1),10);
temp_geneids = nan*ones(length(rnafile)-1,1);
for k=2:length(rnafile),
  if mod(k,1000) == 0,
    k
  end;
  ensgname = rnafile{k}{1};
  ensgid = str2double(ensgname(5:end));
  temp_geneids(k-1) = ensgid;
  % I = find(bininfo(:,5)==ensgid);
  % if length(I) > 0,
  %   for l=1:10,
  %     rna_summaryseries(I,l) = str2double(rnafile{k}{1+l});
  %   end;
  % end;
end;
temp_genenames = cell(length(temp_geneids),1);
for k=1:length(temp_geneids),
  temp_genenames{k} = sprintf('ENSG%011d',temp_geneids(k));
end;


[I, A, B] = intersect(temp_geneids, bininfo(:, 5));
%symbols = r.symbols(A);
symbols = temp_genenames(A);
end;

figdir = '/share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed/newfits';

load /share/synergy/analyses/2012-05-07/pol2rna_joint_fits_2012-05-07.mat

modelok=zeros(length(results_jointmodels),1);
for k=1:length(results_jointmodels),
  if ~isempty(results_jointmodels{k}),
    modelok(k)=1;
  end;
end;
myI = find(modelok==1);

cd(figdir);

for k=1:length(myI),
    l=myI(k);
    plot_model(res.results_jointmodels{l}, res.results_ensemblids(l), res.results_loglikelihoods(l), symbols{l});

    % mkdir(figdir, sprintf('%d', k));
    % cd(sprintf('%s/%d', figdir, k));
    % fprintf('%d/%d\n', find(myI==k), length(myI));
    % res = load(sprintf(filestem, k));
    % I = (res.results_geneindices ~= 0);
    % for l=find(I),
    %    plot_model(res.results_jointmodels{l}, res.results_ensemblids(l), res.results_loglikelihoods(l), symbols{l});
    % end
end
