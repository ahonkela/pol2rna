resultdir = '~/projects/pol2rnaseq/analyses/hmc_results/joint/';
%resultdir = '/share/synergy/analyses/hmc_results/joint/';

d = dir([resultdir '*_samples_2013-09-04_unif0.mat']);

filenames = {};
[filenames{1:length(d),1}] = deal(d.name);
%filenames = filenames(cellfun('length', filenames) == 44);

%N_GOOD = 5;

ids = zeros(size(filenames));
for k=1:length(filenames),
  ids(k) = str2num(filenames{k}(5:15));
end

Isampl = 501:1000;
Isampl_thin = 501:10:1000;

MYP = [2.5 25 50 75 97.5];

if 0,
  load('~/mlprojects/pol2rnaseq/matlab/difficult_gene_files.mat');
  s = struct2cell(genefiles);
  filenames = [filenames; cat(1, s{:})];
  numgenes = length(filenames);
else
  numgenes = length(filenames);
end

%parI = [1:4, 6, 8:10];
parI = [1:6, 8:10];
mygene = '';
genes = cell(length(filenames), 1);
means = zeros(length(filenames), length(parI));
stds = zeros(length(filenames), length(parI));
thetameans = zeros(length(filenames), 1);
prcts = zeros(length(filenames), length(parI), length(MYP));
curI = 1:length(Isampl_thin);
Iincr = length(curI);
h = cell(length(filenames), 1);
for k=1:length(filenames),
  fprintf('%d/%d\n', k, length(filenames));
  r = load([resultdir, filenames{k}]);
  genes{k} = r.gene_name;
  if r.finished,
    hk = cat(3, r.HMCsamples{:});
    hk = hk(Isampl, parI, :);
    hk = permute(hk, [1 3 2]);
    sz = size(hk);
    hk = reshape(hk, [sz(1)*sz(2), sz(3)]);
    prcts(k, :, :) = prctile(hk, MYP)';
    means(k, :) = mean(hk);
    stds(k, :) = std(hk);
  end
end

% p30 = zeros(length(filenames), 1);
% bound30 = sigmoidabTransform(30, 'xtoa', [0, 299]);

% for k=1:length(h),
%   p30(k) = mean(h{k}(:, 5) > bound30);
% end

% baddata = zeros(numgenes, 1);
% Rhat = zeros(numgenes, length(parI));
% thetatruemeans = zeros(numgenes, 1);
% for k=1:numgenes,
%   I = (k-1)*N_GOOD + (1:N_GOOD);
%   S = sum(bsxfun(@minus, means(I,:), median(means(I,:))) .^ 2, 2);
%   J = I(S < 100 & sum(stds(I, :), 2) > 0.1);
%   if length(J) < N_GOOD,
%     baddata(k) = 1;
%     disp(J)
%   end
%   N = length(Isampl);
%   M = length(I);
%   W = mean(stds(I,:).^2);
%   B = N*var(means(I,:));
%   varHatPlus = (N-1)/N * W + 1/N * B;
%   Rhat(k,:) = sqrt(varHatPlus./W);
%   thetatruemeans(k) = mean(thetameans(I));
% end

% K = max(Rhat, [], 2) < 1.2;
% goodgenes = genes(find(K)*5);
% fp = fopen('finished_genes_latest.txt', 'w');
% fprintf(fp, '%s\n', goodgenes{:});
% fclose(fp);

% Kbad = ~(max(Rhat, [], 2) < 1.2);
% badgenes = genes(find(Kbad)*5);
% fp = fopen('unfinished_genes_latest.txt', 'w');
% fprintf(fp, '%s\n', badgenes{:});
% fclose(fp);

% save analysis_dump

% for k=find(max(Rhat, [], 2) > 1.2)',
%   I = (k-1)*5 + (1:5);
%   fprintf('k=%d\n', k);
%   J = Rhat(k,:) > 1.2;
%   disp(Rhat(k,:));
%   disp(means(I,:));
%   disp(stds(I,:));
%   % for k=1:length(I),
%   %   subplot(5, 1, k);
%   %   plot(h{I(k)}(:,J));
%   % end
%   % pause
% end

save result_summary_2013-09-04 means stds prcts genes

pr = load('early_profiles_2013-09-04.mat');
diff1 = max(pr.mu(:, 1:20), [], 2) - min(pr.mu(:, 1:20), [], 2);
diff2 = max(pr.mu(:, 21:30), [], 2) - min(pr.mu(:, 21:30), [], 2);
qrtl = prctile(diff1-diff2, [25, 75]);
bound = qrtl(2) + 1.5 * (qrtl(2) - qrtl(1));
early_dev = (diff1 - diff2 - qrtl(2)) / (qrtl(2) - qrtl(1));

[~, A, B] = intersect(genes, pr.genes);
assert(all(strcmp(genes(A), pr.genes)));
early_devs = zeros(length(genes), 1);
early_devs(A) = early_dev;

delaypcts = sigmoidabTransform(squeeze(prcts(:,5,:)), 'atox', [0, 299]);
fp = fopen('hmc_results_to_browser_2013-09-04.txt', 'w');
fprintf(fp, 'gene 2.5%% 25%% 50%% 75%% 97.5%% begdev\n');
for k=1:length(genes),
  fprintf(fp, '%s %f %f %f %f %f %f\n', genes{k}, delaypcts(k,:), early_devs(k));
end
fclose(fp);
