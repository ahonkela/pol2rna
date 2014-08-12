function plot_sampled_predictions_file(gene, sampledir, filestem, format, SQRTTIME, FILESPEC, PLOTPREMRNA, PREMRNAMODEL, SHORTHIST),

if nargin < 4,
  format = 'png';
end

if nargin < 5,
  SQRTTIME = 1;
end

if nargin < 7,
  PLOTPREMRNA = 0;
end

if nargin < 8,
  PREMRNAMODEL = 0;
end

if nargin < 9,
  SHORTHIST = 0;
end

resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/';
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

d = dir([resultdir sampledir gene FILESPEC]);
filenames = {};
[filenames{1:length(d),1}] = deal(d.name);

Isampl = 501:10:1000;

if length(filenames)~=1,
  fprintf('Found %d != 1 files for gene %s, aborting...\n', ...
          length(filenames), gene);
  return;
end

r = load([resultdir, sampledir, filenames{1}]);

if ~r.finished,
  fprintf('Gene %s not finished yet, exitting...\n', gene);
  return;
end

if PLOTPREMRNA,
  if PREMRNAMODEL,
    load('~/projects/pol2rnaseq/data/pol2_summaryseries_2013_01_02.mat');
    load('~/projects/pol2rnaseq/data/bininfo_dec2012_corrected.mat', 'bininfo');
    pol2I = find(ensg2int(gene) == bininfo(:, 5));
    my_pol2 = pol2_summaryseries(pol2I, :);
  else
    premrna_data = load('~/projects/pol2rnaseq/data/rpkm.mat');
    premrnaI = find(strcmp(gene, premrna_data.genes));
    my_premrna = premrna_data.mT2(premrnaI, :);
  end
end

assert(strcmp(r.gene_name, gene));
mysamples = cat(3, r.HMCsamples{:});
mysamples = mysamples(Isampl, :, :);
mysamples = permute(mysamples, [1 3 2]);
sz = size(mysamples);
mysamples = reshape(mysamples, [sz(1)*sz(2), sz(3)]);

% Check if delay is fixed and only plot delay histogram if it is not
g = gpnddisimLogLikeGradients(r.m);
NODELAYHIST = (g(5) == 0);

if isfield(aliases, gene),
  gene = sprintf('%s (%s)', gene, aliases.(gene));
end

if PREMRNAMODEL,
  if PLOTPREMRNA,
    plot_sampled_predictions_premrna(r.m, mysamples, gene, NODELAYHIST, SQRTTIME, my_pol2, SHORTHIST)
  else
    plot_sampled_predictions_premrna(r.m, mysamples, gene, NODELAYHIST, SQRTTIME, [], SHORTHIST)
  end
else
  if PLOTPREMRNA,
    plot_sampled_predictions(r.m, mysamples, gene, NODELAYHIST, SQRTTIME, my_premrna, SHORTHIST)
  else
    plot_sampled_predictions(r.m, mysamples, gene, NODELAYHIST, SQRTTIME, [], SHORTHIST)
  end
end

set(gcf, 'PaperUnits', 'inches')
switch format,
  case 'png',
    set(gcf, 'PaperPosition', [0, 0, PNG_SIZE/DPI])
    print('-dpng', sprintf('-r%d', DPI), plotfile);
  case 'eps',
    set(gcf, 'PaperPositionMode', 'auto');
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf, 'PaperPosition', [0, 0, 5.0, 5.0])
    print('-depsc2', plotfile);
end
