%function tempanswer=run_all_genetests(startpercent);

tempanswer=1;

timeshift = 300;

if ~exist('initializationtype', 'var'),
    error('No initializationtype set.')
end

if ~exist('dataset', 'var'),
    error('Dataset not selected.')
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

  modelnames = {'joint_nodelay', 'pol2', 'rna'};
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
    fname = sprintf('hmc_results/%s/%s_samples_%s_init%d.mat', ...
		    modelnames{k}, gene_name, id, initializationtype);
    if exist(fname, 'file'),
      fprintf('File %s exists, skipping...\n', fname);
      skipme = 1;
    else
      [m,temptransforminfo]=createNdSimDisim_celltimes_newdata2(...
          timeCell,dataVals,lengthscale,initializationtype,[],[],1, rnaVars);

      oldparams = modelExtractParam(m);
      oldparams(oldparams < -5) = -5;
      oldparams(oldparams > 5) = 5;
      m = modelExpandParam(m, oldparams);

      if 0,
        m.kern.comp{1}.comp{1}.priors{1}.sd = 20;
        m.kern.comp{1}.comp{1}.priors{2}.sd = 20;
        m.kern.comp{1}.comp{2}.priors{1}.sd = 20;
        m.kern.comp{1}.comp{2}.priors{2}.sd = 20;
      end

      m.fix.index = 5;
      m.fix.value = -100;
      
      options = hmcDefaultOptions;
      burninopts = options;
    end
    
    if ~skipme,
      EPS_SCHEDULE = [0.00001 0.0001 0.001 0.003 0.005 0.01 0.03 0.05 0.07 0.1];
      burninopts.epsilon = EPS_SCHEDULE(1);
      goodeps = burninopts.epsilon;
      for myeps = 1:length(EPS_SCHEDULE),
        burninHMCsamples = gpnddisimSampleHMC(m, 0, 100, burninopts);
        if mean(all(diff(burninHMCsamples) == 0, 2)) > 0.2,
          fprintf('Too high rejection rate for gene %s, backtracking eps schedule...\n', gene_name);
          burninopts.epsilon = goodeps;
        else
          goodeps = burninopts.epsilon;
          burninopts.epsilon = EPS_SCHEDULE(myeps);
        end
        oldparams = burninHMCsamples(end,:);
        m = modelExpandParam(m, burninHMCsamples(end,:));
        disp(oldparams)
      end
      fprintf('Burn-in done\n');

      options.epsilon = goodeps;
      % Apply HMC sampling
      HMCsamples = gpnddisimSampleHMC(m, 0, nHMCiters, options);
      HMCsamples = HMCsamples(10:10:end, :);
      save(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples');
    end
  end
end;
