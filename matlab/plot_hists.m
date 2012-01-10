function plot_hists(samples, transforms, names, varargin),

if nargin >= 2,
  for k=1:9,
    samples(:, k) = sigmoidabTransform(samples(:, k), 'atox', transforms{k});
  end
end

for k=1:9,
  subplot(3, 3, k);
  hist(samples(:, k), varargin{:});
  title(names{k})
end
