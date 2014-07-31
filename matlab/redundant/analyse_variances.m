[d, genes] = load_bitseq_data;
[y, yvar] = unlog_bitseq_data(d);
r = load('~/projects/pol2rnaseq/data/rna_new_data2.mat');
r2 = load('~/projects/pol2rnaseq/data/rna_new_data3.mat');
subplot(1, 3, 1);
hist(mean(r.normcounts ./ r.normsds, 2), 50);
title('edgeR');
xlabel('mean / sd')
subplot(1, 3, 2);
hist(mean(y ./ sqrt(yvar), 2), 50);
title('BitSeq');
xlabel('mean / sd')
subplot(1, 3, 3);
hist(mean(r2.normcounts ./ r2.normsds, 2), 50);
title('edgeR new');
xlabel('mean / sd')
