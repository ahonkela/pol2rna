function f = multipriorObjective(params, model)

% MULTIPRIOROBJECTIVE Wrapper function for MULTIPRIOR objective.
% FORMAT
% DESC returns the negative log likelihood of a Gaussian process
% model for single input motifs given the model structure and
% a vector parameters. This allows the use of NETLAB minimisation
% functions to find the model parameters.
% ARG params : the parameters of the model for which the objective
% will be evaluated.
% ARG model : the model structure for which the objective will be
% evaluated.
% RETURN f : the negative log likelihood of the MULTIPRIOR model.
%
% SEEALSO : scg, conjgrad, multipriorCreate, multipriorGradient, multipriorLogLikelihood, multipriorOptimise
% 
% COPYRIGHT : Neil D. Lawrence, 2005, 2006
%
% COPYRIGHT : Jaakko Peltonen, 2011
%
% COPYRIGHT : Antti Honkela, 2007-2012

try,
    model = multipriorExpandParam(model, params);
    f = - multipriorLogLikelihood(model);
catch
    f = -Inf;
end
