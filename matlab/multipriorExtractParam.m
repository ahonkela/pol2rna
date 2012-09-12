function [params, names] = multipriorExtractParam(model);

% MULTIPRIOREXTRACTPARAM Extract params from a multiprior model.
%
% COPYRIGHT : Antti Honkela, 2012

if model.marginalise,
  if nargout < 2,
    params = modelExtractParam(model.basemodel);
  else
    [params, names] = modelExtractParam(model.basemodel);
  end
else
  if nargout < 2,
    params = [model.theta, modelExtractParam(model.basemodel)];
  else,
    [bparams, bnames] = modelExtractParam(model.basemodel);
    params = [model.theta, bparams];
    names{1} = 'Prior selector';
    names(2:model.numParams) = bnames;
  end
end
