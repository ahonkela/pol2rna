timeshift = 300;
mybasedir_data='~/projects/pol2rnaseq/';
mybasedir_analyses=mybasedir_data;
datadir=[mybasedir_data 'data/'];

pol2rnaseqToolboxes;

load([datadir, 'bininfo_nov2014_corrected.mat'], 'bininfo');
load([datadir, 'pol2_summaryseries_2014_11_19.mat']);
normfacts = importdata([datadir, 'rna_norm_factors.txt']);
r = load([datadir, 'info_gene_mean_var.mat']);
act = importdata([datadir, 'BF_RNASeq_gene.txt']);
act_mrna = act.textdata(find(act.data > 3));

pol2act = importdata([datadir, 'Pol2_BF8.txt']);
act_pol2 = pol2act.textdata(find(pol2act.data > 3));

act_genes = unique([act_mrna; act_pol2]);

[I, A, B] = intersect(ensg2int(r.geneID), bininfo(:, 5));

J = zeros(length(act_genes), 1);
for k=1:length(act_genes),
  JJ = strcmp(act_genes{k}, r.geneID(A));
  if any(JJ),
    J(k) = find(JJ);
  end
end

interestinggenes_rna = A(sort(J(find(J))));
interestinggenes_pol2 = B(sort(J(find(J))));

genes = importdata('final_genes.txt');

runids = zeros(size(genes))';
for k=1:length(runids),
  runids(k) = find(interestinggenes_rna == find(strcmp(r.geneID, genes{k})));
end

rna_index = interestinggenes_rna(runids); 
gene_index= interestinggenes_pol2(runids);
gene_name = r.geneID(rna_index);

for k=1:length(runids),
  assert(ensg2int(gene_name{k}) == bininfo(gene_index(k), 5), ...
         'RNA and Pol2 data index mismatch %s != %d', ...
         gene_name{k}, bininfo(gene_index(k), 5));
end

dataVals1=pol2_summaryseries(gene_index,:)';
dataVals2=bsxfun(@rdivide, r.mu(rna_index,:)', normfacts);
rnaVars=bsxfun(@rdivide, r.v(rna_index,:)', (normfacts.^2));
timevector=[0 5 10 20 40 80 160 320 640 1280]' + timeshift;
save ode_data_2015-05-07.mat dataVals1 dataVals2 rnaVars timevector gene_name
