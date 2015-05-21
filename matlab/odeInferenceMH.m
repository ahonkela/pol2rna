%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;
spline_regparam = 0.92; % cross-validation result

if ~exist('seeds', 'var'),
    error('seeds not set.')
end

if ~exist('use_splines', 'var'),
  use_splines = 0;
end

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
genes = importdata('final_genes.txt');

runids = zeros(size(genes))';
for k=1:length(runids),
  runids(k) = find(interestinggenes_rna == find(strcmp(r.geneID, genes{k})));
end

myI = intersect(myI, 1:length(runids));
%myI = 104;
%myI = 1:10;

for i=runids(myI),
  rna_index = interestinggenes_rna(i); 
  gene_index= interestinggenes_pol2(i);
  gene_name = r.geneID{rna_index};

  assert(ensg2int(gene_name) == bininfo(gene_index, 5), ...
         'RNA and Pol2 data index mismatch %s != %d', ...
         gene_name, bininfo(gene_index, 5));
  
  if isnan(gene_index),
    continue;
  end

  dataVals1=pol2_summaryseries(gene_index,:)';
  dataVals2=r.mu(rna_index,:)' ./ normfacts;
  rnaVars=r.v(rna_index,:)' ./ (normfacts.^2);
  rnascale = max(dataVals2) / 10;
  dataVals2 = dataVals2 / rnascale;
  rnaVars = rnaVars / (rnascale.^2);
  dataVals1 = dataVals1 - min(dataVals1);

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

  myindex = find(i==myI);
  for myseedi=1:length(seeds),
    fprintf('Running gene %d/%d, seed %d/%d: %s\n', myindex, length(myI), myseedi, length(seeds), gene_name);

    fname = sprintf('ode_mcmc_results/%s_samples_%s_spl%d_seed%d.mat', ...
                    gene_name, id, use_splines, seeds(myseedi));
    if exist(fname, 'file'),
      fprintf('File %s exists, skipping...\n', fname);
      continue;
    end

    rng(seeds(myseedi));

    if use_splines,
      fprintf('Using spline fit for Pol-II.\n');
      pfit.spline = csaps(log(timevector - min(timevector) + 5), ...
                          dataVals1, spline_regparam);
      pfit.timeshift = 5;
      [samples, accepts] = mhInference(@(params) odeLikelihoodSpline(pfit, dataVals2, rnaVars, timevector, params, 0, 1), 6, 20000, 100);
      save(fname, 'gene_name', 'gene_index', 'samples', 'accepts', 'pfit');
    else
      fprintf('Using linear interpolation for Pol-II.\n');
      [samples, accepts] = mhInference(@(params) odeLikelihoodNoisy(dataVals1, dataVals2, rnaVars, timevector, params, 0, 1), 6, 20000, 100);
      save(fname, 'gene_name', 'gene_index', 'samples', 'accepts');
    end
  end
  % plot
  %odeLikelihood(dataVals1, dataVals2, rnaVars, timevector, params, 1)
  %drawnow;
end;
