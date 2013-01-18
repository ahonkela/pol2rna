function plot_sampled_predictions_file(gene, filestem, format),

if nargin < 3,
  format = 'png';
end

resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
plotdir = '~/projects/pol2rnaseq/analyses/hmc_results/plots/';
aliases = load('~/projects/pol2rnaseq/data/aliases.mat');
aliases = aliases.aliases;

switch format,
  case 'png',
    PNG_SIZE = [640 480];
    DPI = 72;
    plotfile = [plotdir gene filestem '.png'];
  case 'eps',
    plotfile = [plotdir gene filestem '.eps'];
end

if exist(plotfile, 'file'),
  fprintf('Image %s exists, exitting...\n', plotfile);
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

if isfield(aliases, gene),
  gene = sprintf('%s (%s)', gene, aliases.(gene));
end

plot_sampled_predictions(r.m, mysamples, gene)
set(gcf, 'PaperUnits', 'inches')
switch format,
  case 'png',
    set(gcf, 'PaperPosition', [0, 0, PNG_SIZE/DPI])
    print('-dpng', sprintf('-r%d', DPI), plotfile);
  case 'eps',
    set(gcf, 'PaperPosition', [0, 0, 4, 2])
    print('-depsc2', plotfile);
end
