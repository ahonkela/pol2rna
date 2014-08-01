resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
id = '2013-08-30';

d = dir([resultdir '*_samples_' id '_unif0.mat']);

filenames = {};
[filenames{1:length(d),1}] = deal(d.name);

gpids = zeros(size(filenames));
for k=1:length(filenames),
  gpids(k) = str2num(filenames{k}(5:15));
end

load('~/projects/pol2rnaseq/data/bininfo_dec2012_corrected.mat', 'bininfo')
load('~/projects/pol2rnaseq/data/pol2_summaryseries_2013_01_02.mat')       
good = ~(any(isnan(pol2_summaryseries), 2) | sum(pol2_summaryseries == 0, 2) > 1);
pol2ids = sort(bininfo(good,5));
ids = intersect(gpids, pol2ids);
pol2ensg = cell(size(ids));
fp = fopen('results/good_pol2_genes.txt', 'w');
for k=1:length(ids),
  pol2ensg{k} = sprintf('ENSG%011d', ids(k));
  fprintf(fp, '%s\n', pol2ensg{k});
end
fclose(fp);
