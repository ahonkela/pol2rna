%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

if ~exist('initializationtype', 'var'),
    error('No initializationtype set.')
end

if ~exist('do_marginalise', 'var'),
  do_marginalise = 1;
end

%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
%mybasedir_code='~/synergy_data/tempcodebranch/';
%mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
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

addpath(path1,path2,path3,path4,path5,path6,path7,path8)



cd(datadir)
%load h3k4me3_series.mat
%load series_for_matti_ver3.mat
load('pol2_for_matti_ver3.mat', 'bininfo');
load('pol2_summaryseries_2012_09.mat');
r = load('rna_new_data3.mat');
act = importdata('activeGenes_new.txt');
act_genes = act.textdata(find(act.data));
badgenes = importdata('unfinished_genes_2012-11-06.txt');

cd(analysisdir)

[I, A, B] = intersect(r.geneids, bininfo(:, 5));

J = zeros(length(badgenes), 1);
for k=1:length(badgenes),
  JJ = strcmp(badgenes{k}, r.genes(A));
  if any(JJ),
    J(k) = find(JJ);
  end
end

%interestinggenes = A; %find(r.pvals(A) < 0.1);
%interestinggenes = find(r.pvals(A) < 0.9);
%interestinggenes = find(sum(r.counts(A,:), 2) >= 1000);
interestinggenes = A(sort(J(find(J))));

%keyboard;

nHMCburnin = 30000;
nHMCiters = 10000;

for i=myI,
  % Start sampling from the maximum likelihood joint-model fit
  %model = allresults_jointmodels_joint{Ijointrank(i)};

  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  %model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  %model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  %model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
  rna_index = interestinggenes(i); 
  gene_index= B(interestinggenes(i));
  gene_name = r.genes{rna_index};

  if isnan(gene_index),
    continue;
  end

  fprintf('Running gene %d/%d: %s\n', find(i==myI), length(myI), gene_name);

  randn('seed',bininfo(gene_index,5)+13*initializationtype);
  rand('seed',bininfo(gene_index,5)+1234567+13*initializationtype);
  
  dataVals1=pol2_summaryseries(gene_index,:)';
  dataVals2=r.normcounts(rna_index,:)';
  rnaVars=(r.normsds(rna_index,:)').^2;
  timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
    
  lengthscale=10;

  timeCell = {timevector,timevector};
  dataVals = {tempvals1, tempvals2};
  %initializationtype=5;

  if any(isnan(tempvals1)) || any(isnan(tempvals2)),
    warning('NaN data for gene %s, skipping...\n', gene_name);
    continue;
  end

  modelnames = {'joint', 'pol2', 'rna'};
  %for k=1:3,
  % NOTE: run only joint models
  for k=1:1,
    switch k,
     case 1,
      % Create joint model
      timeCell = {timevector,timevector};
      dataVals = {tempvals1, tempvals2};
     case 2,
      % Create Pol2 model
      timeCell = {timevector,[]};
      dataVals = {tempvals1, []};
     case 3,
      % Create RNA model
      timeCell = {[],timevector};
      dataVals = {[], tempvals2};
    end
    
    %initializationtype=5;

    % for l in `cat ~/mlprojects/pol2rnaseq/matlab/unfinished_genes_latest.txt` ; do for k in ${l}_*.mat ; do cp -p $k ${k/.mat/_orig.mat} ; done ; done
    
    skipme = 0;
    fname = sprintf('hmc_results/%s/%s_samples_%s_init%d.mat', ...
		    modelnames{k}, gene_name, id, initializationtype);
    fname_orig = sprintf('hmc_results/%s/%s_samples_%s_init%d_orig.mat', ...
                         modelnames{k}, gene_name, id, initializationtype);

    origres = load(fname_orig);
    %if exist(fname, 'file'),
    %  fprintf('File %s exists, skipping...\n', fname);
    %  skipme = 1;
    %else
      %try,
      %[m,temptransforminfo]=createNdSimDisim_celltimes_newdata2(...
      %	    timeCell,dataVals,lengthscale,initializationtype,[],[],1, rnaVars);
        %catch,
	%warning('Unable to create %s model for gene %s\n', modelnames{k}, gene_name);
	%skipme = 1;
        %end
        %end

    % Create the multiprior model
    m = origres.m;
    priorspec = 'kern.comp{1}.comp{2}.priors{3}';
    eval(['p = m.', priorspec, ';']);
    priors = {p, p};
    priors{2}.mu = -4;
    priors{2}.sd = 2;
    if do_marginalise,
      m2 = multipriorCreate(m, priorspec, priors, struct('marginalise', 1));
    else
      m2 = multipriorCreate(m, priorspec, priors, struct());
    end
    options = hmcDefaultOptions;
    options.epsilon=0.05;

    m2 = modelExpandParam(m2, origres.HMCsamples(end,:));
    %if ~skipme,
      % Apply HMC sampling
    burninopts = options;
    oldparams = origres.HMCsamples(end,:);
    for myeps = 1:5,
      burninopts.epsilon = 0.01 * myeps;
      burninHMCsamples = multipriorSampleHMC(m2, 0, 100, burninopts);
      if all(oldparams == burninHMCsamples(end,:)),
        fprintf('Unable to get gene %s going...\n', gene_name);
      end
      oldparams = burninHMCsamples(end,:);
      m2 = modelExpandParam(m2, burninHMCsamples(end,:));
    end
    fprintf('Burn-in done\n');
    HMCsamples = multipriorSampleHMC(m2, 0, nHMCiters, options);
    HMCsamples = HMCsamples(10:10:end, :);
    save(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples');
    %end
  end
end;
