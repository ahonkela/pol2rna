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


load('simulated_data.mat');

datasize = size(rnadata);

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
  gene_name = sprintf('Synthetic %d', i);

  if isnan(gene_index),
    continue;
  end

  fprintf('Running gene %d/%d: %s\n', find(i==myI), length(myI), gene_name);

  randn('seed',bininfo(gene_index,5)+13*initializationtype);
  rand('seed',bininfo(gene_index,5)+1234567+13*initializationtype);
  
  dataVals1=pol2data;
  dataVals2=squeeze(rnadata(rna_index(1), rna_index(2), :));
  rnaVars=0.01 * ones(size(dataVals2));
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

    skipme = 0;
    fname = sprintf('hmc_results/%s/%s_samples_%s_init%d.mat', ...
		    modelnames{k}, gene_name, id, initializationtype);
    if exist(fname, 'file'),
      fprintf('File %s exists, skipping...\n', fname);
      skipme = 1;
    else
      [m,temptransforminfo]=createNdSimDisim_celltimes_newdata2(...
          timeCell,dataVals,lengthscale,initializationtype,[],[],1, rnaVars);

      % Hack: do not change the transformed decay even though the range changes
      oldparams = modelExtractParam(m);

      % Tweak RNA decay and Pol2 noise prior
      %trsets = gpnddisimExtractParamTransformSettings(m);
      % Set RNA decay bounds and prior
      %trsets{3} = [5e-4, 0.35];
      %m.kern.comp{1}.comp{2}.priors{1}.mu = -2;
      % Set Pol2 noise bounds
      %lb = min(min(m.kern.comp{2}.comp{2}.fixedvariance), 0.01*var(tempvals1));
      %ub = min(25*max(m.kern.comp{2}.comp{2}.fixedvariance), ...
      %         0.25*var(tempvals1));
      %trsets{6} = [lb, ub];

      % Set lengthscale bounds
      %trsets{1} = 1 ./ [1280, 20].^2;
      %m = gpnddisimExpandParamTransformSettings(m, trsets);

      oldparams(oldparams < -5) = -5;
      oldparams(oldparams > 5) = 5;
      m = modelExpandParam(m, oldparams);

      % Create the multiprior model
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
      oldparams = modelExtractParam(m2);
      options = hmcDefaultOptions;
      burninopts = options;
    end
    
    if ~skipme,
      EPS_SCHEDULE = [0.00001 0.0001 0.001 0.003 0.005 0.01 0.03 0.05];
      burninopts.epsilon = EPS_SCHEDULE(1);
      goodeps = burninopts.epsilon;
      for myeps = 1:length(EPS_SCHEDULE),
        burninHMCsamples = multipriorSampleHMC(m2, 0, 100, burninopts);
        if mean(all(diff(burninHMCsamples) == 0, 2)) > 0.2,
          fprintf('Too high rejection rate for gene %s, backtracking eps schedule...\n', gene_name);
          burninopts.epsilon = goodeps;
        else
          goodeps = burninopts.epsilon;
          burninopts.epsilon = EPS_SCHEDULE(myeps);
        end
        oldparams = burninHMCsamples(end,:);
        m2 = modelExpandParam(m2, burninHMCsamples(end,:));
        disp(oldparams)
      end
      fprintf('Burn-in done\n');

      options.epsilon = goodeps;
      % Apply HMC sampling
      HMCsamples = multipriorSampleHMC(m2, 0, nHMCiters, options);
      HMCsamples = HMCsamples(10:10:end, :);
      save(fname, 'gene_name', 'gene_index', 'm', 'HMCsamples');
    end
  end
end;
