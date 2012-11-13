function save_sampled_predictions_file(gene),

resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
savedir = '~/projects/pol2rnaseq/analyses/hmc_results/profiles/';

savefile = [savedir gene '_hmc_2012-10-08.mat'];

if exist(savefile, 'file'),
  fprintf('Image %s exists, exitting...\n', savefile);
  return;
end

d = dir([resultdir gene '*.mat']);
filenames = {};
[filenames{1:length(d),1}] = deal(d.name);
% exclude init3, as it seems very unreliable
I = cellfun('isempty', strfind(filenames, 'init3'));
filenames = filenames(I);

Isampl = 501:10:1000;
mysamples = zeros(5*length(Isampl), 10);

for k=1:length(filenames),
  r = load([resultdir, filenames{k}]);
  assert(strcmp(r.gene_name, gene));
  mysamples((1:length(Isampl)) + (k-1)*length(Isampl), :) = ...
      r.HMCsamples(Isampl, :);
end

t_pred = (((0:100)/100*sqrt(1280)).^2 + 300)';

r = gpnddisimSamplePredictions(r.m, mysamples, t_pred, 500);
p = prctile(r, [2.5 97.5]);
mu = mean(r);

save(savefile, 'p', 'mu', 'gene');
