function model = multipriorParamInit(model)

% MULTIPRIORPARAMINIT Initialise the parameters of an MULTIPRIOR model.
% FORMAT
% DESC initialises the prior selector to pick the first prior.
% ARG model : the input model to initialise.
% RETURN model : the initialised model.
%
% SEEALSO : modelParamInit, multipriorCreate
%
% COPYRIGHT : Antti Honkela, 2012

if ~model.marginalise,
  model.theta = 1;
end
