function g = multipriorLogLikeGradients(model)

% MULTIPRIORLOGLIKEGRADIENTS Multiprior model gradients.
% FORMAT
% DESC computes the gradients of the log likelihood of a
% multiprior model with respect to the parameters.
% ARG model : the model structure for computing the log likelihood.
% RETURN g : the gradients of the model log likelihood.
%
% SEEALSO : modelLogLikeihood, multipriorLogLikelihood
%
% COPYRIGHT : Antti Honkela, 2012

if ~model.marginalise,
  error('multipriorLogLikeGradients unsupported for non-marginalised models');
else
  g = zeros(1, model.numParams);
  for k=1:model.numPriors,
    eval(['model.basemodel.', model.priorarg, '=model.priors{k};']);
    g = g + modelLogLikeGradients(model.basemodel);
  end
end
