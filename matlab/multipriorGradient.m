function g = multipriorGradient(params, model)

% GPDISIMGRADIENT Gradient wrapper for a MULTIPRIOR model.
% FORMAT 
% DESC wraps the log likelihood gradient function to return the
% gradient of the negative of the log likelihood. This can then be
% used in, for example, NETLAB, minimisation tools.
% ARG params : the parameters of the model.
% ARG model : the model for which gradients will be computed.
% RETURN g : the returned gradient of the negative log likelihood
% for the given parameters.
%
% SEEALSO : scg, conjgrad, multipriorCreate, multipriorObjective, multipriorLogLikeGradient
%
% COPYRIGHT : Neil D. Lawrence, 2006
%
% COPYRIGHT : Antti Honkela, 2007-2012

try,
    model = multipriorExpandParam(model, params);
    g = - multipriorLogLikeGradients(model);
catch,
    g = 0*multipriorLogLikeGradients(model);
end
