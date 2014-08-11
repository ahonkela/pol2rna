function save_sampled_predictions_file(gene, sampledir, filestem, FILESPEC, SQRTTIME),

resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/';
savedir = '~/projects/pol2rnaseq/analyses/hmc_results/profiles/';

savefile = [savedir gene filestem '.mat'];

if exist(savefile, 'file'),
  fprintf('Image %s exists, exitting...\n', savefile);
  return;
end

fprintf('Finding result: %s\n', [resultdir sampledir gene FILESPEC]);

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

assert(strcmp(r.gene_name, gene));
mysamples = cat(3, r.HMCsamples{:});
mysamples = mysamples(Isampl, :, :);
mysamples = permute(mysamples, [1 3 2]);
sz = size(mysamples);
mysamples = reshape(mysamples, [sz(1)*sz(2), sz(3)]);

switch SQRTTIME,
  case 0,
    t_pred0 = [linspace(0, 157.5, 64), exp(linspace(log(160), log(1280), 61))]';
    t_pred = [(200:5:295)'; (t_pred0 + 300)];
  case 1,
    t_pred = [(200:5:295)'; (((0:100)/100*sqrt(1280)).^2 + 300)'];
  otherwise,
    error('Unsupported SQRTTIME')
end

r = gpnddisimSamplePredictions(r.m, mysamples, t_pred, 500);
p = prctile(r, [2.5 97.5]);
mu = mean(r);

save(savefile, 'p', 'mu', 'gene', 't_pred');
