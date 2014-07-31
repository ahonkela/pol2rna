function dic = compute_dic_for_model(model),

N = size(model.HMCsamples, 1);
m = model.m;
lls = zeros(N, 1);
for k=1:N,
  if rem(k, 100) == 0,
    fprintf('Running iteration %d/%d\n', k, N);
  end
  m = modelExpandParam(m, model.HMCsamples(k, :));
  lls(k) = gpnddisimLogLikelihood(m, 1);
end

Dbar = -2 * mean(lls);
trans = gpnddisimExtractParamTransformSettings(m);
mu = transform_params(mean(transform_params(model.HMCsamples, trans)), trans, 'xtoa');
m = modelExpandParam(m, mu);
Dhat = -2 * gpnddisimLogLikelihood(m, 1);
dic = 2*Dbar - Dhat;



function p = transform_params(p0, trans, type),

if nargin < 3,
  type = 'atox';
end

p = zeros(size(p0));
for k=1:size(p, 2),
  p(:, k) = sigmoidabTransform(p0(:, k), type, trans{k});
end
