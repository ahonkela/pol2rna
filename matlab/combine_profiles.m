resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/profiles/';

%id = '2013-08-30';
id = '2013-11-05';

d = dir([resultdir 'ENSG*' id '.mat']);

filenames = {};
[filenames{1:length(d),1}] = deal(d.name);

cwd = pwd();
cd(resultdir);
r = load(filenames{1});
genes = cell(length(filenames), 1);
mu = zeros(length(filenames), length(r.mu));
prctiles = zeros(length(filenames), length(r.mu), 2);
t_pred = r.t_pred;
for k=1:length(filenames),
  if rem(k, 100)==0,
    fprintf('%d/%d\n', k, length(filenames));
  end
  r = load(filenames{k});
  genes{k} = r.gene;
  mu(k,:) = r.mu;
  prctiles(k,:,:) = r.p';
end
save(['all_profiles_' id '.mat'], 'genes', 'mu', 'prctiles', 't_pred');
cd(cwd);
