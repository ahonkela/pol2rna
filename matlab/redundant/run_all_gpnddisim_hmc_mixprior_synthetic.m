%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

if ~exist('initializationtype', 'var'),
    error('No initializationtype set.')
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

if ~exist('use_uniform_priors', 'var'),
  use_uniform_priors = 1;
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


load('simulated_data_2013-08-09.mat');

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

  randn('seed',gene_index+13*initializationtype);
  rand('seed',gene_index+1234567+13*initializationtype);
  
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
  %initializationtype=5;

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
    
    %initializationtype=5;

    skipme = 0;
    fname = sprintf('hmc_results/%s/%s_samples_%s_unif%d_init%d.mat', ...
		    modelnames{k}, gene_name, id, use_uniform_priors, ...
                    initializationtype);
    if exist(fname, 'file'),
      fprintf('File %s exists, skipping...\n', fname);
      skipme = 1;
    else
      [m,temptransforminfo]=createNdSimDisim_celltimes_newdata3(...
          timeCell,dataVals,lengthscale,initializationtype,[],[],1, rnaVars,1,use_uniform_priors);

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

      m.fix(1).index = 10;
      if use_uniform_priors,
        m.fix(1).value = 1e-100;
      else
        m.fix(1).value = -100;
      end

      if nodelay,
        m.fix(2).index = 5;
        if use_uniform_priors,
          m.fix(2).value = 1e-100;
        else
          m.fix(2).value = -100;
        end
      end

      options = hmcDefaultOptions;
    end
    
    if ~skipme,
      [options.epsilon, options.scales, m] = gpnddisimTuneHMC(m, options);
      fprintf('Adaptation done\n');

      % Apply HMC sampling
      HMCsamples = gpnddisimSampleHMC(m, 1, nHMCiters, options);
      HMCsamples = HMCsamples(10:10:end, :);
      save(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples', 'options');
    end
  end
end;
