function compute_r2_file(gene, sampledir, filestem),

resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/';
savedir = '~/projects/pol2rnaseq/analyses/hmc_results/r2_results/';

savefile = [savedir gene filestem '.mat'];

if exist(savefile, 'file'),
  fprintf('Image %s exists, exitting...\n', savefile);
  return;
end

d = dir([resultdir sampledir gene '*.mat']);
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

ss_res = zeros(size(mysamples, 1), 1);
ss_tot = zeros(size(mysamples, 1), 1);
for k=1:size(mysamples, 1),
  m = modelExpandParam(r.m, mysamples(k,:,:));
  ss_res(k) = sum(((m.K(11:20,1:10) / m.K(1:10,1:10) * m.m(1:10)) - m.m(11:20)).^2);
  ss_tot(k) = sum((m.m(11:20) - m.mu(11:20)).^2);
end
%ss_tot = sum((m.m(11:20) - mean(m.m(11:20))).^2);
r2 = 1 - mean(ss_res ./ ss_tot);

save(savefile, 'r2', 'ss_res', 'ss_tot', 'gene');
