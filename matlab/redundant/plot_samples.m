function plot_samples(samples, thin),

if nargin < 2,
  thin = 1;
end

samples = samples(thin:thin:end, :, :);
d = size(samples);

s = permute(samples, [2 1 3]);
s = reshape(s, [d(2), d(1)*d(3)]);
plot(s');

%legend(1:d(2))
