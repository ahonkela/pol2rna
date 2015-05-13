genes = {'ENSG00000162949',
         'ENSG00000019549',
         'ENSG00000163659',
         'ENSG00000117298',
         'ENSG00000166483',
         'ENSG00000181026',
         'ENSG00000140961'};
seeds = 11:14;
id = '2015-05-06';

res = cell(length(genes), length(seeds));
means = zeros(length(genes), length(seeds), 5);
medians = zeros(length(genes), length(seeds), 5);
stds = zeros(length(genes), length(seeds), 5);
for k=1:length(genes),
  for l=1:length(seeds),
    fname = sprintf('ode_mcmc_results/%s_samples_%s_seed%d.mat', genes{k}, id, seeds(l));
    res{k,l} = load(fname);
    means(k, l, :) = mean(res{k,l}.samples(101:end, :));
    medians(k, l, :) = median(res{k,l}.samples(101:end, :));
    stds(k, l, :) = std(res{k,l}.samples(101:end, :));
  end
end

