%id = '2013-08-30';
%id = '2013-11-05';
id = 'final';
profiledir='~/projects/pol2rnaseq/analyses/hmc_results/profiles/';
d = dir([profiledir 'ENSG*_' id '.mat']);
filenames = {};
[filenames{1:length(d),1}] = deal(d.name);

mu = zeros(length(filenames), 30);
genes = cell(size(filenames));
for k=1:length(filenames),
  if rem(k, 100)==0,
    fprintf('%d/%d\n', k, length(filenames));
  end
  rr = load([profiledir filenames{k}]);
  mu(k,:) = rr.mu(1:30);
  genes{k} = filenames{k}(1:15);
end
save(['results/early_profiles_' id '.mat'], 'mu', 'genes');
