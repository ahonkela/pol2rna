function model = multipriorExpandParam(model, params);

% MULTIPRIOREXPANDPARAM Update multiprior model with vector of parameters.
%
% COPYRIGHT : Antti Honkela, 2012

if model.marginalise,
  model.basemodel = modelExpandParam(model.basemodel, params);
else
  model.theta = params(1);
  model.basemodel = modelExpandParam(model.basemodel, params(2:end));
  eval(['model.basemodel.', model.priorarg, '=model.priors{model.theta};']);
end
