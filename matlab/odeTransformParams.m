function boundedparams = odeTransformParams(params),

% beta0
% beta^2
% alpha
% Delta
% m0
% sigma_n
bounds = [0, 1;
          1e-6, 1;
          1e-6, log(2);
          0, 299;
          0, 10;
          0, 10];

assert(size(params, 2) == 5 || size(params, 2) == 6);
boundedparams = zeros(size(params));
for k=1:size(params, 2),
  boundedparams(:, k) = bounds(k, 1) + (bounds(k, 2)-bounds(k, 1))* sigmoid(params(:, k));
end
