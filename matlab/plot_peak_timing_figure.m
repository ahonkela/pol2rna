% addpath ~/mlprojects/matlab/general
% importLatest('pol2rnaseq')
% pol2rnaseqToolboxes

rna = load('~/projects/pol2rnaseq/data/rpkm.mat');
pol2 = load('~/projects/pol2rnaseq/data/pol2_summaryseries_2013_01_02.mat');
load('~/projects/pol2rnaseq/data/bininfo_dec2012_corrected.mat');
pol2.bininfo = bininfo;

rna_act = importdata('~/projects/pol2rnaseq/data/BF_RNASeq_gene.txt');

assert(all(strcmp(rna.genes, rna_act.textdata)));

rna.bfs = rna_act.data;

rna.ids = ensg2int(rna.genes);
[~, rnamax] = max(rna.mT2, [], 2);
[~, mrnamax] = max(rna.mT1, [], 2);
[~, polmax] = max(pol2.pol2_summaryseries, [], 2);
[I, A, B] = intersect(rna.ids, pol2.bininfo(:, 5));
J = rna.bfs(A) > 3;
indices = [polmax(B(J)), rnamax(A(J))];
peaks = accumarray(indices, 1, [10 10]);

indices2 = [polmax(B(J)), mrnamax(A(J))];
peaks2 = accumarray(indices2, 1, [10 10]);

sum(diag(peaks)) / sum(sum(peaks))
(sum(diag(peaks)) + sum(diag(peaks, 1)) + sum(diag(peaks, -1))) / sum(sum(peaks))

exppeaks = sum(peaks, 2) * sum(peaks, 1) / sum(sum(peaks));

sums = zeros(1, 19);
expsums = zeros(1, 19);
for k=-9:9,
    sums(k+10) = sum(diag(peaks, k));
    expsums(k+10) = sum(diag(exppeaks, k));
end

exppeaks2 = sum(peaks2, 2) * sum(peaks2, 1) / sum(sum(peaks2));

sums2 = zeros(1, 19);
expsums2 = zeros(1, 19);
for k=-9:9,
    sums2(k+10) = sum(diag(peaks2, k));
    expsums2(k+10) = sum(diag(exppeaks2, k));
end

plot(-9:9, sums / sum(sums), 'k');
hold on
plot(-9:9, sums2 / sum(sums2), 'k--');
plot(-9:9, expsums / sum(expsums), 'k:');
hold off
xlabel('Difference in peak timing')
ylabel('Frequency')
legend('Intronic-Pol2', 'mRNA-Pol2', 'Expected (intronic-Pol2)', ...
       'location', 'NorthEast');

set(gcf, 'PaperPosition', [0 0 6 4]);
