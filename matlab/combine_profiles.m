resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/profiles/';

d = dir([resultdir 'ENSG*.mat']);

filenames = {};
[filenames{1:length(d),1}] = deal(d.name);

cd(resultdir);
genes = cell(length(filenames), 1);
mu = zeros(length(filenames), 202);
prctiles = zeros(length(filenames), 202, 2);
for k=1:length(filenames),
  r = load(filenames{k});
  genes{k} = r.gene;
  mu(k,:) = r.mu;
  prctiles(k,:,:) = r.p';
end
save('all_profiles_2013-01-18.mat', 'genes', 'mu', 'prctiles');
