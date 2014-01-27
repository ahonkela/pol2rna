function samples = multipriorSampleHMC(model, display, iters, options);

% MULTIPRIORSAMPLEHMC Do HMC sampling for the MULTIPRIOR model.
% FORMAT
% DESC performs HMC sampling for the Gaussian process single input
% motif model for a given number of iterations.
% ARG model : the model to be optimised.
% ARG display : whether or not to display while optimisation
% proceeds, set to 2 for the most verbose and 0 for the least
% verbose.
% ARG iters : number of samples to return.
% RETURN samples : the samples.
%
% SEEALSO : hmc, gpsimCreate, gpsimGradient, gpsimObjective
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2012
%
% COPYRIGHT : Jaakko Peltonen, 2011


if nargin < 3
  iters = 2000;
  if nargin < 2
    display = 1;
  end
end


[params,paramnames] = modelExtractParam(model);

if nargin < 4,
  options = hmcDefaultOptions;
end
if display
  options.verbose = display;
end

samples = zeros(iters, length(params));

if model.marginalise,
  hmcopt = optOptions;
  hmcopt(1) = options.verbose;
  if length(params) <= 100
    hmcopt(9) = 1;
  end

  % Momentum persistence
  hmcopt(5) = 1;

  % Leapfrog steps
  hmcopt(7) = options.tau;

  % hmcopt(18) = epsilon, step length
  hmcopt(18) = options.epsilon;

  % Number of samples to return
  hmcopt(14) = iters;
  
  samples = hmc('modelObjective', params, hmcopt, ...
                'modelGradient', model);
else
  for k=1:iters,
    try,
      params(2:end) = hmcStep('modelObjective', params(2:end),  options, ...
                              'modelGradient', model.basemodel);
    catch me,
      warning(['Error in hmcStep: ', me.message]);
    end
    model = modelExpandParam(model, params);
    params(1) = multipriorSampleTheta(model, options);
    samples(k,:) = params;
    eval(['model.basemodel.', model.priorarg, '=model.priors{params(1)};']);
  end
end



function theta = multipriorSampleTheta(model, options),

lls = zeros(1, model.numPriors);
for k=1:model.numPriors,
  eval(['model.basemodel.', model.priorarg, '=model.priors{k};']);
  lls(k) = modelLogLikelihood(model.basemodel);
end
if options.verbose,
  disp(lls)
end
lls = lls + log(model.thetaprior);
lls = lls - logsumexp(lls);
% Shortcut to sample from discrete distribution without Statistics Toolbox
[~,theta] = histc(rand(1),[0;cumsum(exp(lls(:)))]);

eval(['model.basemodel.', model.priorarg, '=model.priors{theta};']);
