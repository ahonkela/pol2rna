datadir='~/projects/pol2rnaseq/data/';

times = [0 5 10 20 40 80 160 320 640 1280];

% pol-II
load([datadir 'bininfo_nov2014_corrected.mat'], 'bininfo');
load([datadir 'pol2_summaryseries_2014_11_19.mat']);

[~, I_pol2] = sort(bininfo(:, 5));

for k=1:length(times),
  fprintf('pol-II, k=%d\n', k);
  fp = fopen(sprintf('data/MCF7_PolII_t%04d.txt', times(k)), 'w');
  fprintf(fp, 'gene_id\tpolII_density\n');
  for l0=1:size(bininfo, 1),
    l = I_pol2(l0);
    if ~any(isnan(pol2_summaryseries(l,:))),
      fprintf(fp, 'ENSG%011d\t%f\n', bininfo(l, 5), pol2_summaryseries(l, k));
    end
  end
  fclose(fp);
end

% RNA
normfacts = importdata([datadir 'rna_norm_factors.txt']);
r = load([datadir 'info_gene_mean_var.mat']);

[~, I_mrna] = sort(r.geneID);

for k=1:length(times),
  fprintf('mRNA, k=%d\n', k);
  fp = fopen(sprintf('data/MCF7_mRNA_t%04d.txt', times(k)), 'w');
  fprintf(fp, 'gene_id\tgene_name\texpression_mean\texpression_sd\n');
  for l0=1:length(r.geneID),
    l = I_mrna(l0);
    fprintf(fp, '%s\t%s\t%f\t%f\n', r.geneID{l}, r.gene_name{l}, r.mu(l, k)/normfacts(k), sqrt(r.v(l, k))/normfacts(k));
  end
  fclose(fp);
end
