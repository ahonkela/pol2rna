id = '2013-08-30';

profiles = load(['~/projects/pol2rnaseq/analyses/hmc_results/profiles/all_profiles_' id '.mat']);

%t_pred = (((0:100)/100*sqrt(1280)).^2)';
T_MIN = 300;
t_pred = profiles.t_pred - T_MIN;
realt = find(t_pred >= 0);

mygenes = importdata('../R/analysed_genes.txt');

[I, A, B] = intersect(mygenes, profiles.genes);

dataVals1 = profiles.mu(B, realt);
dataVals2 = profiles.mu(B, realt + length(profiles.t_pred));
dataVals1(dataVals1(:) < 0) = 0;

sz = size(dataVals1);
dataVals1 = dataVals1 ./ repmat(max(dataVals1, [], 2), [1, sz(2)]);
dataVals2 = dataVals2 ./ repmat(max(dataVals2, [], 2), [1, sz(2)]);

val1mean = mean(dataVals1 .* repmat(1:sz(2), [sz(1), 1]), 2);
val2mean = mean(dataVals2 .* repmat(1:sz(2), [sz(1), 1]), 2);

[~, I] = sort(val1mean + val2mean);
subplot(1, 2, 1);
imagesc(dataVals1(I, :))
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(I, :))
title('mRNA')
print -depsc2 -opengl profiles_heatmap_meansort

[~, I] = sort(val1mean);
subplot(1, 2, 1);
imagesc(dataVals1(I, :))
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(I, :))
title('mRNA')
print -depsc2 -opengl profiles_heatmap_pol2sort

[~, I] = sort(val2mean);
subplot(1, 2, 1);
imagesc(dataVals1(I, :))
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(I, :))
title('mRNA')
print -depsc2 -opengl profiles_heatmap_mrnasort

Y = pdist([dataVals1, dataVals2]);
Z = linkage(Y, 'average');
order = optimalleaforder(Z, Y);

subplot(1, 2, 1);
imagesc(dataVals1(order, :))
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(order, :))
title('mRNA')
print -depsc2 -opengl profiles_heatmap_hclustsort
