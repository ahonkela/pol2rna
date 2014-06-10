resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
%resultdir = '/share/synergy/analyses/hmc_results/joint/';

d = dir([resultdir '*.mat']);

filenames = {};
[filenames{1:length(d),1}] = deal(d.name);

N_GOOD = 5;

ids = zeros(size(filenames));
for k=1:length(filenames),
  ids(k) = str2num(filenames{k}(5:15));
end
goodids = find(accumarray(ids, 1) == N_GOOD);
Ifiles = ismember(ids, goodids);

Isampl = 501:1000;
Isampl_thin = 501:10:1000;

MYP = [2.5 25 50 75 9.75];

filenames = filenames(Ifiles);

mygene = '';
genes = cell(length(filenames), 1);
means = zeros(length(filenames), 9);
stds = zeros(length(filenames), 9);
thetameans = zeros(length(filenames), 1);
prcts = zeros(length(filenames)/N_GOOD, 9, length(MYP));
hgene = zeros(length(Isampl)/2, 9);
curI = 1:length(Isampl_thin);
Iincr = length(curI);
%h = cell(length(filenames), 1);
for k=1:length(filenames),
  fprintf('%d/%d\n', k, length(filenames));
  r = load([resultdir, filenames{k}]);
  genes{k} = r.gene_name;
  hk = r.HMCsamples(Isampl, [1:6, 8:10]);
  pp = [normpdf(hk(:, 5), 0, 2), normpdf(hk(:, 5), -4, 2)];
  thetameans(k) = sum(mean(pp ./ repmat(sum(pp, 2), [1, 2])) .* [1 2]);
  if ~strcmp(r.gene_name, mygene),
    if k>1,
      prcts(floor(k/N_GOOD), :, :) = prctile(hgene, MYP)';
    end
    curI = 1:length(Isampl_thin);
    mygene = r.gene_name;
  end
  hgene(curI,:) = r.HMCsamples(Isampl_thin, [1:6, 8:10]);
  curI = curI + Iincr;
  means(k, :) = mean(hk);
  stds(k, :) = std(hk);
end
prcts(end, :, :) = prctile(hgene, MYP)';

% p30 = zeros(length(filenames), 1);
% bound30 = sigmoidabTransform(30, 'xtoa', [0, 299]);

% for k=1:length(h),
%   p30(k) = mean(h{k}(:, 5) > bound30);
% end

baddata = zeros(length(goodids), 1);
Rhat = zeros(length(goodids), 9);
for k=1:length(goodids),
  I = (k-1)*N_GOOD + (1:N_GOOD);
  S = sum(bsxfun(@minus, means(I,:), median(means(I,:))) .^ 2, 2);
  I = I(S < 100 & sum(stds(I, :), 2) > 0.1);
  if length(I) < N_GOOD,
    baddata(k) = 1;
    disp(I)
  end
  N = length(Isampl);
  M = length(I);
  W = mean(stds(I,:).^2);
  B = N*var(means(I,:));
  varHatPlus = (N-1)/N * W + 1/N * B;
  Rhat(k,:) = sqrt(varHatPlus./W);
end