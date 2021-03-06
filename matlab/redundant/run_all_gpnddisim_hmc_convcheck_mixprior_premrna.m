%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

if ~exist('seeds', 'var'),
    error('seeds not set.')
end

if ~exist('nodelay', 'var'),
  nodelay = 0;
end

if ~exist('use_pol2_fixedvar', 'var'),
  use_pol2_fixedvar = 0;
end

if ~exist('use_uniform_priors', 'var'),
  use_uniform_priors = 0;
end

if ~exist('do_not_normalise', 'var'),
  do_not_normalise = 0;
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
%load('pol2_for_matti_ver3.mat', 'bininfo');
%load('pol2_summaryseries_2012_09.mat');
%r = load('rna_new_data4.mat');
%act = importdata('activeGenes_new.txt');
normfacts = importdata('rna_norm_factors.txt');
r = load('info_gene_mean_var.mat');

premrna_data = load('rpkm.mat');

% Assume both data sets have the same genes
assert(all(strcmp(r.geneID, premrna_data.genes)));

act = importdata('BF_RNASeq_gene.txt');
act_mrna = act.textdata(find(act.data > 3));

pol2act = importdata('Pol2_BF8.txt');
act_pol2 = pol2act.textdata(find(pol2act.data > 3));

act_genes = unique([act_mrna; act_pol2]);

%extra_genes = {'ENSG00000064195'; 'ENSG00000187109'; 'ENSG00000198668'; 'ENSG00000142178'; 'ENSG00000183484'; 'ENSG00000105372'; 'ENSG00000168140'; 'ENSG00000134107'; 'ENSG00000166226'; 'ENSG00000101331'; 'ENSG00000164626'; 'ENSG00000099860'; 'ENSG00000010278'; 'ENSG00000137193'; 'ENSG00000181026'; 'ENSG00000163659'};

%act_genes = sort([act_mrna; extra_genes]);

if exist('run_just_bad', 'var'),
  act_genes = importdata(run_just_bad);
end

cd(analysisdir)

[~, ~, J] = intersect(act_genes, r.geneID);
assert(length(J) == length(intersect(r.geneID(J), act_genes)));
interestinggenes_rna = J;

%interestinggenes = A; %find(r.pvals(A) < 0.1);
%interestinggenes = find(r.pvals(A) < 0.9);
%interestinggenes = find(sum(r.counts(A,:), 2) >= 1000);

nHMCiters = 10000;

myI = intersect(myI, 1:length(interestinggenes_rna));

for i=myI,
  % Start sampling from the maximum likelihood joint-model fit
  %model = allresults_jointmodels_joint{Ijointrank(i)};

  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  %model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  %model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  %model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
  rna_index = interestinggenes_rna(i); 
  gene_index= rna_index;
  gene_name = r.geneID{rna_index};

  fprintf('Running gene %d/%d: %s\n', find(i==myI), length(myI), gene_name);

  dataVals1=premrna_data.mT2(gene_index,:)' ./ normfacts;
  premrnaVars = premrna_data.vT2(gene_index,:)' ./ (normfacts.^2);
  dataVals2=r.mu(rna_index,:)' ./ normfacts;
  rnaVars=r.v(rna_index,:)' ./ (normfacts.^2);
  timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;

  temptimes=timevector;
  tempvals1=dataVals1;
  tempvals2=dataVals2;
    
  lengthscale=10;

  timeCell = {timevector,timevector};
  dataVals = {tempvals1, tempvals2};
  dataVars = {premrnaVars, rnaVars};

  if any(isnan(tempvals1)) || any(isnan(tempvals2)),
    warning('NaN data for gene %s, skipping...\n', gene_name);
    continue;
  end

  if nodelay,
    modelnames = {'joint_nodelay', 'pol2', 'rna'};
  else
    modelnames = {'joint', 'pol2', 'rna'};
  end
  %for k=1:3,
  % NOTE: run only joint models
  for k=1:1,
    switch k,
     case 1,
      % Create joint model
      timeCell = {timevector,timevector};
      dataVals = {tempvals1, tempvals2};
      dataVars = {premrnaVars, rnaVars};
     case 2,
      % Create Pol2 model
      timeCell = {timevector,[]};
      dataVals = {tempvals1, []};
      dataVars = {premrnaVars, []};
     case 3,
      % Create RNA model
      timeCell = {[],timevector};
      dataVals = {[], tempvals2};
      dataVars = {[], rnaVars};
    end
    
    fname = sprintf('hmc_results/%s/%s_samples_%s_unif%d.mat', ...
		    modelnames{k}, gene_name, id, use_uniform_priors);
    if exist(fname, 'file'),
      fprintf('File %s exists, loading existing results...\n', fname);
      load(fname);
    else
      
      [m,temptransforminfo]=createNdSimDisim_celltimes_newdata3(...
          timeCell,dataVals,lengthscale,11,[],[],1, dataVars,do_not_normalise,use_uniform_priors);

      oldparams = modelExtractParam(m);
      %oldparams(oldparams < -5) = -5;
      %oldparams(oldparams > 5) = 5;
      if use_uniform_priors,
	bounds = gpnddisimExtractParamTransformSettings(m);
	bounds = cat(1, bounds{:})';
	I1 = (oldparams < bounds(1,:));
	oldparams(I1) = bounds(1, I1);
	I2 = (oldparams > bounds(2,:));
	oldparams(I2) = bounds(2, I2);

	m = modelExpandParam(m, oldparams);
      end

      if use_uniform_priors,
        m.kern.comp{1}.comp{2}.priors{3} = ...
            priorCreate('mixture', struct('types', {{'truncatedGamma', 'uniform'}}));
        m.kern.comp{1}.comp{2}.priors{3}.index = 5;
        m.kern.comp{1}.comp{2}.priors{3} = ...
            priorSetBounds(m.kern.comp{1}.comp{2}.priors{3}, [0, 299]);
        m.kern.comp{1}.comp{2}.priors{3}.comp{1} = ...
            priorExpandParam(m.kern.comp{1}.comp{2}.priors{3}.comp{1}, [0, log(0.1)]);
      else
        m.kern.comp{1}.comp{2}.priors{3} = ...
            priorCreate('mixture', struct('types', {{'logisticNormal', 'logisticNormal'}}));
        m.kern.comp{1}.comp{2}.priors{3}.index = 5;
        m.kern.comp{1}.comp{2}.priors{3} = ...
            priorSetBounds(m.kern.comp{1}.comp{2}.priors{3}, [0, 299]);
        m.kern.comp{1}.comp{2}.priors{3}.comp{1} = ...
            priorExpandParam(m.kern.comp{1}.comp{2}.priors{3}.comp{1}, [-2, log(2)]);
        m.kern.comp{1}.comp{2}.priors{3}.comp{2} = ...
            priorExpandParam(m.kern.comp{1}.comp{2}.priors{3}.comp{2}, [-2, log(2)]);
      end

      %m.fix(1).index = 10;
      %if use_uniform_priors,
      %  m.fix(1).value = 1e-100;
      %else
      %  m.fix(1).value = -100;
      %end

      if nodelay,
        m.fix(1).index = 5;
        if use_uniform_priors,
          m.fix(1).value = 1e-100;
        else
          m.fix(1).value = -100;
        end
      end

      options = hmcDefaultOptions;
      options.maxit = 100000;
      finished = 0;
      samples_done = 0;
      states.rng = cell(size(seeds));
      HMCsamples = cell(size(seeds));
      states.params = cell(size(seeds));
      tuning.epsilons = zeros(size(seeds));
      tuning.scales = cell(size(seeds));
      for chain=1:length(seeds),
        rng(gene_index+13*seeds(chain));
      
	if use_uniform_priors,
          % Random initial parameters
          initparams = bounds(1,:) + ...
		       rand(size(oldparams)) .* (bounds(2,:)-bounds(1,:));

          % Initialise to a small delay
          initparams(5) = min(initparams(5), 10*rand(1));
	else
          % Random initial parameters
	  initparams = randn(size(oldparams));

          % Initialise to a small delay
          initparams(5) = initparams(5) - 2;
	end

        m = modelExpandParam(m, initparams);
      
        [tuning.epsilons(chain), tuning.scales{chain}, m] = gpnddisimTuneHMC(m, options);
        fprintf('Adaptation done\n');

        states.params{chain} = modelExtractParam(m);
        % Save chain random state
        states.rng{chain} = rng;
      end
      [~,bestI] = min(tuning.epsilons);
      options.epsilon = tuning.epsilons(bestI(1));
      options.scales = tuning.scales{bestI(1)};
      save(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples', 'options', ...
           'finished', 'states', 'samples_done', 'tuning');
    end     % if no previous save file found

    while ~finished && samples_done < options.maxit,
      for chain=1:length(seeds),
        % Restore chain random state
        rng(states.rng{chain});
        m = modelExpandParam(m, states.params{chain});
        
        % Apply HMC sampling
        HMCsamples{chain} = gpnddisimSampleHMC(m, 0, nHMCiters, options);
        states.rng{chain} = rng;
        states.params{chain} = HMCsamples{chain}(end, :);
        HMCsamples{chain} = HMCsamples{chain}(10:10:end, :);
      end

      g = gpnddisimLogLikeGradients(m);
      Irealparams = find(g);
      samples_done = samples_done + nHMCiters;
      N = 500;
      means = cellfun(@(x) squeeze(mean(x(end-N:end,:), 1)), HMCsamples, 'UniformOutput', false);
      stds = cellfun(@(x) squeeze(std(x(end-N:end,:), 0, 1)), HMCsamples, 'UniformOutput', false);
      W = mean(cat(1, stds{:}).^2);
      B = N*var(cat(1, means{:}));
      varHatPlus = (N-1)/N * W + 1/N * B;
      Rhat = sqrt(varHatPlus./W);
      fprintf('%d samples done, Rhat:\n', samples_done);
      disp(Rhat);
      finished = all(Rhat(Irealparams) < 1.2);
      
      safeSave(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples', 'options', ...
               'finished', 'states', 'samples_done', 'tuning');
    end
  end
end;
