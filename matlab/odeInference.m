%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

% if ~exist('seeds', 'var'),
%     error('seeds not set.')
% end

% if ~exist('nodelay', 'var'),
%   nodelay = 0;
% end

% if ~exist('use_pol2_fixedvar', 'var'),
%   use_pol2_fixedvar = 0;
% end

% if ~exist('use_uniform_priors', 'var'),
%   use_uniform_priors = 0;
% end

% if ~exist('do_not_normalise', 'var'),
%   do_not_normalise = 0;
% end

mybasedir_data='~/projects/pol2rnaseq/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
datadir=[mybasedir_data 'data/'];
analysisdir=[mybasedir_analyses 'analyses/'];

pol2rnaseqToolboxes;

addpath /home/fs/ahonkela/src/bayesopt/matlab

if ~exist('bininfo', 'var'),
  load([datadir, 'bininfo_nov2014_corrected.mat'], 'bininfo');
  load([datadir, 'pol2_summaryseries_2014_11_19.mat']);
  normfacts = importdata([datadir, 'rna_norm_factors.txt']);
  r = load([datadir, 'info_gene_mean_var.mat']);
  act = importdata([datadir, 'BF_RNASeq_gene.txt']);
  act_mrna = act.textdata(find(act.data > 3));

  pol2act = importdata([datadir, 'Pol2_BF8.txt']);
  act_pol2 = pol2act.textdata(find(pol2act.data > 3));

  act_genes = unique([act_mrna; act_pol2]);
end

if exist('run_just_bad', 'var'),
  act_genes = importdata(run_just_bad);
end

%cd(analysisdir)

[I, A, B] = intersect(ensg2int(r.geneID), bininfo(:, 5));

J = zeros(length(act_genes), 1);
for k=1:length(act_genes),
  JJ = strcmp(act_genes{k}, r.geneID(A));
  if any(JJ),
    J(k) = find(JJ);
  end
end

%interestinggenes = A; %find(r.pvals(A) < 0.1);
%interestinggenes = find(r.pvals(A) < 0.9);
%interestinggenes = find(sum(r.counts(A,:), 2) >= 1000);
interestinggenes_rna = A(sort(J(find(J))));
interestinggenes_pol2 = B(sort(J(find(J))));

nHMCiters = 10000;

genes = {'ENSG00000162949',
         'ENSG00000019549',
         'ENSG00000163659',
         'ENSG00000117298',
         'ENSG00000166483',
         'ENSG00000181026',
         'ENSG00000140961'};
myI = zeros(size(genes))';
for k=1:length(myI),
  myI(k) = find(interestinggenes_rna == find(strcmp(r.geneID, genes{k})));
end

%myI = intersect(myI, 1:length(interestinggenes_rna));
%myI = 104;
%myI = 1:10;

fitparams = cell(size(myI));
lls = zeros(size(myI));

for i=myI,
  % Start sampling from the maximum likelihood joint-model fit
  %model = allresults_jointmodels_joint{Ijointrank(i)};

  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  %model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  %model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  %model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
  rna_index = interestinggenes_rna(i); 
  gene_index= interestinggenes_pol2(i);
  gene_name = r.geneID{rna_index};

  assert(ensg2int(gene_name) == bininfo(gene_index, 5), ...
         'RNA and Pol2 data index mismatch %s != %d', ...
         gene_name, bininfo(gene_index, 5));
  
  if isnan(gene_index),
    continue;
  end

  myindex = find(i==myI);
  fprintf('Running gene %d/%d: %s\n', myindex, length(myI), gene_name);

  dataVals1=pol2_summaryseries(gene_index,:)';
  dataVals2=r.mu(rna_index,:)' ./ normfacts;
  rnaVars=r.v(rna_index,:)' ./ (normfacts.^2);
  timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
    
  lengthscale=10;

  timeCell = {timevector,timevector};
  dataVals = {tempvals1, tempvals2};

  if any(isnan(tempvals1)) || any(isnan(tempvals2)),
    warning('NaN data for gene %s, skipping...\n', gene_name);
    continue;
  end

  dataVals1 = dataVals1 - min(dataVals1);
  %  params0 = [-3, -3, -3, 0, dataVals2(1)];
  %params0 = [-10, -2.5, -3, 3.5, 1];
  %options = optimset('MaxFunEvals', 10000, 'MaxIter', 10000, 'display', 'iter');

  boparams.n_iterations = 190;
  boparams.n_init_samples = 10;
  boparams.crit_name = 'cEI';
  boparams.surr_name = 'sStudentTProcessNIG';
  boparams.noise = 1e-6;
  boparams.kernel_name = 'kMaternARD5';
  boparams.kernel_hp_mean = [1];
  boparams.kernel_hp_std = [10];
  boparams.verbose_level = 1;
  boparams.log_filename = 'matbopt.log';
  
  % beta0, beta, alpha, Delta, m(0)
  %lb = [0, 0, 1, 0, 0];
  %ub = [1, 5, 300, 120, max(dataVals2)];
  lb = [-5, -5, -5, -8, -5];
  ub = [5, 5, 5, 0, 5];
  %[params, ll] = fminsearch(@(params) -odeLikelihood(dataVals1, dataVals2, rnaVars, timevector, params), params0, options);
  [params, ll] = bayesoptcont(@(params) -odeLikelihood(dataVals1, dataVals2, rnaVars, timevector, params), 5, boparams, lb, ub);
  fitparams{myindex} = params;
  lls(myindex) = ll;
  
  % plot
  odeLikelihood(dataVals1, dataVals2, rnaVars, timevector, params, 1)
  drawnow;
end;
