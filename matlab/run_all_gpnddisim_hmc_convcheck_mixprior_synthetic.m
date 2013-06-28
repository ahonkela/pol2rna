%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

if ~exist('seeds', 'var'),
    error('seeds not set.')
end

if ~exist('dataset', 'var'),
    error('Dataset not selected.')
end

if ~exist('nodelay', 'var'),
  nodelay = 0;
end

if ~exist('use_pol2_fixedvar', 'var'),
  use_pol2_fixedvar = 0;
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


load('simulated_data.mat');

datasize = size(rnadata{dataset});

nHMCiters = 10000;

myI = intersect(myI, 1:prod(datasize(1:2)));

for i=myI,
  % Start sampling from the maximum likelihood joint-model fit
  %model = allresults_jointmodels_joint{Ijointrank(i)};

  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  %model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  %model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  %model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
  [rna_index(1), rna_index(2)] = ind2sub(datasize(1:2), i); 
  gene_index = 1;
  gene_name = sprintf('Synthetic_%02d', i);

  if isnan(gene_index),
    continue;
  end

  fprintf('Running gene %d/%d: %s\n', find(i==myI), length(myI), gene_name);

  if ndims(pol2data{dataset}) > 2,
    dataVals1=squeeze(pol2data{dataset}(rna_index(1), rna_index(2), :));
  else
    dataVals1=pol2data{dataset};
  end
  dataVals2=squeeze(rnadata{dataset}(rna_index(1), rna_index(2), :));
  if use_pol2_fixedvar,
    rnaVars = {0.01 * ones(size(dataVals1)), 0.01 * ones(size(dataVals2))};
  else
    rnaVars=0.01 * ones(size(dataVals2));
  end
  timevector=t_gen{dataset}' + timeshift;

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
     case 2,
      % Create Pol2 model
      timeCell = {timevector,[]};
      dataVals = {tempvals1, []};
     case 3,
      % Create RNA model
      timeCell = {[],timevector};
      dataVals = {[], tempvals2};
    end
    
    fname = sprintf('hmc_results/%s/%s_samples_%s.mat', ...
		    modelnames{k}, gene_name, id);
    if exist(fname, 'file'),
      fprintf('File %s exists, loading existing results...\n', fname);
      load(fname);
    else
      
      [m,temptransforminfo]=createNdSimDisim_celltimes_newdata3(...
          timeCell,dataVals,lengthscale,11,[],[],1, rnaVars,1,1);

      oldparams = modelExtractParam(m);
      %oldparams(oldparams < -5) = -5;
      %oldparams(oldparams > 5) = 5;
      bounds = gpnddisimExtractParamTransformSettings(m);
      bounds = cat(1, bounds{:})';
      I1 = (oldparams < bounds(1,:));
      oldparams(I1) = bounds(1, I1);
      I2 = (oldparams > bounds(2,:));
      oldparams(I2) = bounds(2, I2);

      m = modelExpandParam(m, oldparams);

      m.kern.comp{1}.comp{2}.priors{3} = ...
          priorCreate('mixture', struct('types', {{'truncatedGamma', 'uniform'}}));
      m.kern.comp{1}.comp{2}.priors{3}.index = 5;
      m.kern.comp{1}.comp{2}.priors{3} = ...
          priorSetBounds(m.kern.comp{1}.comp{2}.priors{3}, [0, 299]);
      m.kern.comp{1}.comp{2}.priors{3}.comp{1} = ...
          priorExpandParam(m.kern.comp{1}.comp{2}.priors{3}.comp{1}, [0, log(0.1)]);

      if nodelay,
        m.fix.index = 5;
        m.fix.value = 1e-100;
      end

      g = gpnddisimLogLikeGradients(m);
      Irealparams = find(g);

      options = hmcDefaultOptions;
      options.maxit = 100000;
      finished = 0;
      samples_done = 0;
      states.rand = cell(size(seeds));
      states.randn = cell(size(seeds));
      HMCsamples = cell(size(seeds));
      states.params = cell(size(seeds));
      tuning.epsilons = zeros(size(seeds));
      tuning.scales = cell(size(seeds));
    
      for chain=1:length(seeds),
        randn('seed',gene_index+13*seeds(chain));
        rand('seed',gene_index+1234567+13*seeds(chain));
      
        % Random initial parameters
        initparams = bounds(1,:) + ...
            rand(size(oldparams)) .* (bounds(2,:)-bounds(1,:));

        % Initialise to a small delay
        initparams(5) = min(initparams(5), 10*rand(1));
        m = modelExpandParam(m, initparams);
      
        [tuning.epsilons(chain), tuning.scales{chain}, m] = gpnddisimTuneHMC(m, options);
        fprintf('Adaptation done\n');

        states.params{chain} = modelExtractParam(m);
        % Save chain random state
        states.rand{chain} = rand('state');
        states.randn{chain} = randn('state');
      end
      [~,bestI] = max(tuning.epsilons);
      options.epsilon = tuning.epsilons(bestI(1));
      options.scales = tuning.scales{bestI(1)};
      save(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples', 'options', ...
           'finished', 'states', 'samples_done', 'tuning');
    end     % if no previous save file found

    while ~finished && samples_done < options.maxit,
      for chain=1:length(seeds),
        % Restore chain random state
        rand('state', states.rand{chain});
        randn('state', states.randn{chain});
        m = modelExpandParam(m, states.params{chain});
        
        % Apply HMC sampling
        HMCsamples{chain} = gpnddisimSampleHMC(m, 0, nHMCiters, options);
        states.rand{chain} = rand('state');
        states.randn{chain} = randn('state');
        states.params{chain} = modelExtractParam(m);
        HMCsamples{chain} = HMCsamples{chain}(10:10:end, :);
      end

      samples_done = samples_done + nHMCiters;
      means = cellfun(@(x) squeeze(mean(x(end-500:end,:), 1)), HMCsamples, 'UniformOutput', false);
      stds = cellfun(@(x) squeeze(std(x(end-500:end,:), 0, 1)), HMCsamples, 'UniformOutput', false);
      W = cellfun(@(x) mean(x.^2, 2), stds, 'UniformOutput', false);
      B = cellfun(@(x) N*var(x, 0, 2), means, 'UniformOutput', false);
      varHatPlus = cellfun(@(W, B) (N-1)/N * W + 1/N * B, W, B, 'UniformOutput', false);
      Rhat = cellfun(@(a, b) sqrt(a ./ b), varHatPlus, W);
      fprintf('%d samples done, Rhat:\n', samples_done);
      disp(Rhat);
      finished = all(Rhat(Irealparams) < 1.2);
      
      safeSave(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples', 'options', ...
               'finished', 'states', 'samples_done', 'tuning');
    end
  end
end;
