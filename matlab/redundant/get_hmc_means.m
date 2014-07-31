% analyse_hmc

parI = [1:4, 6, 8:10];

[~, paramnames] =  modelExtractParam(r.m);
trsets = gpnddisimExtractParamTransformSettings(r.m);
trsets = trsets(parI);
paramnames = paramnames(parI);
genemeans = squeeze(mean(reshape(means, [5, size(means, 1)/5, length(parI)])));
truemeans = zeros(size(genemeans));
for k=1:length(trsets),
  truemeans(:, k) = sigmoidabTransform(genemeans(:, k), 'atox', trsets{k});
end

[~, J] = sort(goodgenes);

load('~/projects/pol2rnaseq/data/aliases.mat');

profiles = load('~/projects/pol2rnaseq/analyses/hmc_results/profiles/all_profiles_2013-01-18.mat');

t_pred = (((0:100)/100*sqrt(1280)).^2)';
[~, polmax] = max(profiles.mu(:, 1:101), [], 2);
poltmax = t_pred(polmax);

corrs = load('~/projects/pol2rnaseq/data/pol2intronicrna_correlation.mat');
I = corrs.pol2intronicrna_correlation > 0.5;
corrs.pol2intronicrna_correlation = corrs.pol2intronicrna_correlation(I,:);
corrs.bininfo = corrs.bininfo(I,:);
corrs.genes = cell(size(corrs.pol2intronicrna_correlation));
for k=1:length(corrs.genes),
  corrs.genes{k} = sprintf('ENSG%011d', corrs.bininfo(k,5));
end
[~, A, B] = intersect(corrs.genes, profiles.genes);

mygenes = corrs.genes(A);
pol2intcorr = corrs.pol2intronicrna_correlation(A);
poltmax = poltmax(B);

assert(all(strcmp(mygenes, profiles.genes(B))));

[~, A2, B2] = intersect(goodgenes, profiles.genes(B));
assert(all(B2' == 1:length(B)));
assert(all(strcmp(mygenes, goodgenes(A2))));

delaymeans = truemeans(A2, 5);

fid = fopen('pol2max_and_delays.txt', 'w');
fprintf(fid, 'ENSG\tgene\ttmax\tdelay\n');
for k=1:length(mygenes),
  if isfield(aliases, mygenes{k}),
    fprintf(fid, '%s\t%s\t%f\t%f\n', mygenes{k}, aliases.(mygenes{k}), poltmax(k), delaymeans(k));
  else
    %fprintf(fid, '%s\tN/A\t%f\t%f\n', mygenes{k}, poltmax(k), delaymeans(k));
  end
end
fclose(fid);
