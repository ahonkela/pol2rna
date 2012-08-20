function plot_sampled_predictions_file(gene, format),

if nargin < 2,
  format = 'png';
end

resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
plotdir = '~/projects/pol2rnaseq/analyses/hmc_results/plots/';

switch format,
  case 'png',
    PNG_SIZE = [640 320];
    DPI = 72;
    plotfile = [plotdir gene '_hmc_2012-08-03.png'];
  case 'eps',
    plotfile = [plotdir gene '_hmc_2012-08-03.eps'];
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
