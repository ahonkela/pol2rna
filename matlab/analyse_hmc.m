
d = dir('~/projects/pol2rnaseq/analyses/hmc_results/*.mat');

genes = cell(length(d), 1);
means = zeros(length(d), 9);
stds = zeros(length(d), 9);
h = cell(length(d), 1);
for k=1:length(d),
  fprintf('%d/%d\n', k, length(d));
  r = load(['~/projects/pol2rnaseq/analyses/hmc_results/', d(k).name]);
  genes{k} = r.gene_name;
  h{k} = r.HMCsamples(501:end, [1:6, 8:10]);
  means(k, :) = mean(h{k});
  stds(k, :) = std(h{k});
end

p30 = zeros(length(d), 1);
bound30 = sigmoidabTransform(30, 'xtoa', [0, 299]);

for k=1:length(h),
  p30(k) = mean(h{k}(:, 5) > bound30);
end
