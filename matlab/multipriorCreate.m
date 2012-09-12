function model = multipriorCreate(basemodel, priorarg, priors, options)

% MULTIPRIORCREATE Create a multiple prior model.

model.type = 'multiprior';
model.basemodel = basemodel;
model.priorarg = priorarg;
model.priors = priors;
model.numPriors = length(priors);
model.thetaprior = 1/model.numPriors * ones(1, model.numPriors);
if isfield(options, 'marginalise') && options.marginalise,
  model.marginalise = 1;
  model.numParams = basemodel.numParams;
else
  model.marginalise = 0;
  model.numParams = basemodel.numParams + 1;
end

model = multipriorParamInit(model);
