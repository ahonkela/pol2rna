function ll = multipriorLogLikelihood(model)

% MULTIPRIORLOGLIKELIHOOD Multiprior model log likelihood.
% FORMAT
% DESC computes the log likelihood of a multiprior model. 
% ARG model : the model structure for computing the log likelihood.
% RETURN ll : the model log likelihood.
%
% SEEALSO : modelLogLikeihood
%
% COPYRIGHT : Antti Honkela, 2012

if ~model.marginalise,
  error('multipriorLogLikelihood unsupported for non-marginalised models');
else
  lls = zeros(1, model.numPriors);
  for k=1:model.numPriors,
    eval(['model.basemodel.', model.priorarg, '=model.priors{k};']);
    lls(k) = modelLogLikelihood(model.basemodel);
  end
  lls = lls + log(model.thetaprior);
  ll = sum(lls);
end
