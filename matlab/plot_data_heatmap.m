datadir = '~/projects/pol2rnaseq/data/';

%load h3k4me3_series.mat
%load series_for_matti_ver3.mat
%load('pol2_for_matti_ver3.mat', 'bininfo');
%load('pol2_summaryseries_2012_09.mat');
load([datadir, 'bininfo_dec2012_corrected.mat'], 'bininfo');
load([datadir 'pol2_summaryseries_2013_01_02.mat']);
%r = load('rna_new_data4.mat');
%act = importdata('activeGenes_new.txt');
normfacts = importdata([datadir 'rna_norm_factors.txt']);
r = load([datadir 'info_gene_mean_var.mat']);

mygenes = importdata('../R/analysed_genes.txt');

[I, A, B] = intersect(ensg2int(r.geneID), bininfo(:, 5));

J = zeros(length(mygenes), 1);
for k=1:length(mygenes),
  JJ = strcmp(mygenes{k}, r.geneID(A));
  if any(JJ),
    J(k) = find(JJ);
  end
end

%interestinggenes = A; %find(r.pvals(A) < 0.1);
%interestinggenes = find(r.pvals(A) < 0.9);
%interestinggenes = find(sum(r.counts(A,:), 2) >= 1000);
interestinggenes_rna = A(sort(J(find(J))));
interestinggenes_pol2 = B(sort(J(find(J))));

dataVals1=pol2_summaryseries(interestinggenes_pol2,:);
dataVals2=r.mu(interestinggenes_rna,:) ./ repmat(normfacts', [length(interestinggenes_rna), 1]);

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
print -depsc2 -opengl data_heatmap_meansort

[~, I] = sort(val1mean);
subplot(1, 2, 1);
imagesc(dataVals1(I, :))
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(I, :))
title('mRNA')
print -depsc2 -opengl data_heatmap_pol2sort

[~, I] = sort(val2mean);
subplot(1, 2, 1);
imagesc(dataVals1(I, :))
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(I, :))
title('mRNA')
print -depsc2 -opengl data_heatmap_mrnasort

Y = pdist([dataVals1, dataVals2]);
Z = linkage(Y, 'average');
order = optimalleaforder(Z, Y);

subplot(1, 2, 1);
imagesc(dataVals1(order, :))
title('Pol-II')
subplot(1, 2, 2);
imagesc(dataVals2(order, :))
title('mRNA')
print -depsc2 -opengl data_heatmap_hclustsort
