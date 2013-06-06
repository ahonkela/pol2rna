normfacts = importdata('~/projects/pol2rnaseq/data/rna_norm_factors.txt');
r = load('~/projects/pol2rnaseq/data/info_gene_mean_var.mat');
act = importdata('~/projects/pol2rnaseq/data/BF_RNASeq_gene.txt');
assert(all(strcmp(act.rowheaders, r.geneID)))
rnadata = r.mu(act.data > 3,:) ./ repmat(normfacts', [sum(act.data > 3), 1]);
[~, maxdrop] = min(diff(rnadata, [], 2), [], 2);

t = [0 5 10 20 40 80 160 320 640 1280];

t_half = repmat(diff(t), [size(rnadata, 1), 1]) .* log(2) ./ log(rnadata(:, 2:end) ./ rnadata(:, 1:end-1));
t_half(t_half < 0) = Inf;
t_half_est = min(t_half(:, 1:end), [], 2);
hist(log(t_half_est(~isinf(t_half_est))) / log(10), 100)
ylabel('# of genes')
xlabel('upper bound of log_{10}(t_{1/2})')
set(gcf, 'PaperPosition', [0 0 6 4]);
