% analyse_hmc_new
id = '2013-08-30';

%parI = [1:4, 6, 8:10];

[~, paramnames] =  modelExtractParam(r.m);
trsets = gpnddisimExtractParamTransformSettings(r.m);
trsets = trsets(parI);
paramnames = paramnames(parI);
genemedians = medians;
truemedians = zeros(size(genemedians));
for k=1:length(trsets),
  truemedians(:, k) = sigmoidabTransform(genemedians(:, k), 'atox', trsets{k});
end

load('~/projects/pol2rnaseq/data/aliases.mat');

profiles = load(['~/projects/pol2rnaseq/analyses/hmc_results/profiles/all_profiles_' id '.mat']);

%t_pred = (((0:100)/100*sqrt(1280)).^2)';
T_MIN = 300;
t_pred = profiles.t_pred - T_MIN;
realt = find(t_pred >= 0);
[~, polmax] = max(profiles.mu(:, realt), [], 2);
poltmax = t_pred(realt(polmax));

diff1 = max(profiles.mu(:, 1:20), [], 2) - min(profiles.mu(:, 1:20), [], 2);
diff2 = max(profiles.mu(:, 21:30), [], 2) - min(profiles.mu(:, 21:30), [], 2);
qrtl = prctile(diff1-diff2, [25, 75]);
bound = qrtl(2) + 1.5 * (qrtl(2) - qrtl(1));
begdev10 = (diff1 - diff2 - qrtl(2)) / (qrtl(2) - qrtl(1));

diff1 = max(profiles.mu(:, 1:20), [], 2) - min(profiles.mu(:, 1:20), [], 2);
diff2 = max(profiles.mu(:, 21:36), [], 2) - min(profiles.mu(:, 21:36), [], 2);
qrtl = prctile(diff1-diff2, [25, 75]);
bound = qrtl(2) + 1.5 * (qrtl(2) - qrtl(1));
begdev30 = (diff1 - diff2 - qrtl(2)) / (qrtl(2) - qrtl(1));

corrs = load('~/projects/pol2rnaseq/data/pol2intronicrna_correlation.mat');
%I = corrs.pol2intronicrna_correlation > 0.5;
%corrs.pol2intronicrna_correlation = corrs.pol2intronicrna_correlation(I,:);
%corrs.bininfo = corrs.bininfo(I,:);
corrs.genes = cell(size(corrs.pol2intronicrna_correlation));
for k=1:length(corrs.genes),
  corrs.genes{k} = sprintf('ENSG%011d', corrs.bininfo(k,5));
end
[~, A, B] = intersect(corrs.genes, profiles.genes);

mygenes = corrs.genes(A);
pol2intcorr = corrs.pol2intronicrna_correlation(A);
poltmax = poltmax(B);
begdev10 = begdev10(B);
begdev30 = begdev30(B);

assert(all(strcmp(mygenes, profiles.genes(B))));

goodI = ~all(medians==0, 2);
goodgenes = genes(goodI);
goodmedians = truemedians(goodI, :);

[~, A2, B2] = intersect(goodgenes, profiles.genes(B));
assert(all(B2' == 1:length(B)));
assert(all(strcmp(mygenes, goodgenes(A2))));

delaymedians = goodmedians(A2, 5);

fid = fopen(['pol2max_and_meddelays_' id '.txt'], 'w');
fprintf(fid, 'ENSG\tgene\ttmax\tmeddelay\tcorr\tbegdev10\tbegdev30\n');
for k=1:length(mygenes),
  if ~isfinite(pol2intcorr(k)),
    continue;
  end
  if isfield(aliases, mygenes{k}),
    fprintf(fid, '%s\t%s\t%f\t%f\t%f\t%f\t%f\n', mygenes{k}, aliases.(mygenes{k}), poltmax(k), delaymedians(k), pol2intcorr(k), begdev10(k), begdev30(k));
  else
    fprintf(fid, '%s\tNA\t%f\t%f\t%f\t%f\t%f\n', mygenes{k}, poltmax(k), delaymedians(k), pol2intcorr(k), begdev10(k), begdev30(k));
  end
end
fclose(fid);
