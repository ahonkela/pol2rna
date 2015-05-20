%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

if ~exist('seeds', 'var'),
    error('seeds not set.')
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
genes = { 'ENSG00000163634',
'ENSG00000115641',
'ENSG00000197451',
'ENSG00000185090',
'ENSG00000126756',
'ENSG00000171067',
'ENSG00000166166',
'ENSG00000160058',
'ENSG00000179913',
'ENSG00000197375',
'ENSG00000214530',
'ENSG00000197170',
'ENSG00000107815',
'ENSG00000115875',
'ENSG00000132661',
'ENSG00000110092',
'ENSG00000187145',
'ENSG00000196139',
'ENSG00000173890',
'ENSG00000183578',
'ENSG00000064195',
'ENSG00000181610',
'ENSG00000170442',
'ENSG00000169255',
'ENSG00000126524',
'ENSG00000258289',
'ENSG00000167173',
'ENSG00000131143',
'ENSG00000092010',
'ENSG00000137038',
'ENSG00000120533',
'ENSG00000181026',
'ENSG00000103257',
'ENSG00000131051',
'ENSG00000062725',
'ENSG00000167767',
'ENSG00000072210',
'ENSG00000128567',
'ENSG00000168140',
'ENSG00000198910',
'ENSG00000099992',
'ENSG00000056736',
'ENSG00000125691',
'ENSG00000142864',
'ENSG00000164125',
'ENSG00000048162',
'ENSG00000023445',
'ENSG00000123066',
'ENSG00000173227',
'ENSG00000179388',
'ENSG00000136603',
'ENSG00000164684',
'ENSG00000143126',
'ENSG00000173821',
'ENSG00000164938',
'ENSG00000109452',
'ENSG00000182022',
'ENSG00000165272',
'ENSG00000121671',
'ENSG00000154767',
'ENSG00000103061',
'ENSG00000197622',
'ENSG00000112701',
'ENSG00000120688',
'ENSG00000166949',
'ENSG00000119285',
'ENSG00000254635',
'ENSG00000143553',
'ENSG00000092969',
'ENSG00000127084',
'ENSG00000006652',
'ENSG00000099622',
'ENSG00000144228',
'ENSG00000156127',
'ENSG00000135052',
'ENSG00000248099',
'ENSG00000073282',
'ENSG00000166483',
'ENSG00000071242',
'ENSG00000113971',
'ENSG00000119318',
'ENSG00000096968',
'ENSG00000240498',
'ENSG00000106003',
'ENSG00000162949',
'ENSG00000124831',
'ENSG00000177119',
'ENSG00000168209',
'ENSG00000258890',
'ENSG00000140548',
'ENSG00000163993',
'ENSG00000138119',
'ENSG00000173848',
'ENSG00000135245',
'ENSG00000103035',
'ENSG00000109321',
'ENSG00000184012',
'ENSG00000157483',
'ENSG00000166888',
'ENSG00000143379' };
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

    fname = sprintf('ode_mcmc_results/%s_samples_%s_seed%d.mat', ...
                    gene_name, id, seeds(myseedi));
    if exist(fname, 'file'),
      fprintf('File %s exists, skipping...\n', fname);
      continue;
    end

    rng(seeds(myseedi));
    [samples, accepts] = mhInference(@(params) odeLikelihoodNoisy(dataVals1, dataVals2, rnaVars, timevector, params, 0, 1), 6, 20000, 100);
    save(fname, 'gene_name', 'gene_index', 'samples', 'accepts');
  end
  % plot
  %odeLikelihood(dataVals1, dataVals2, rnaVars, timevector, params, 1)
  %drawnow;
end;
