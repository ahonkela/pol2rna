function [samples, accepts] = mhInference(logp, dim, N, thin),

x = randn(1, dim);
x(4) = -5;

propstd = 0.3;

samples = zeros(N/thin, dim);
accepts = 0;

ll = logp(x);

for k=1:N,
  xp = x + propstd * randn(1, dim);
  lp = logp(xp);
  accp = exp(lp - ll);
  if rand(1) < accp,
    accepts = accepts + 1;
    x = xp;
    ll = lp;
  end
  if mod(k, thin)==0,
    fprintf('%d samples done, acceptance rate %.2f, current logprob %.2f, current state:\n', k, accepts/k, ll);
    disp(x);
    samples(k/thin, :) = x;
  end
end
fprintf('Done, acceptance rate %.2f.\n', accepts/k);
