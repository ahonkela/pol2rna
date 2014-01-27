rpkm=load('/triton/ics/project/synergy/data/Hande_codes/rpkm1.mat')
commondisp = importdata('~/Dropbox/projects/synergy/RNA_analysis/common_dispersion.txt')
normcounts = rpkm.MT1 ./ rpkm.ef_l .* repmat(mean(rpkm.ef_l, 2), [1, 10]);
normsds = sqrt(effcounts + commondisp * effcounts.^2);
genes = rpkm.genes;
geneids=zeros(size(rpkm.genes));
for k=1:length(geneids),
  geneids(k)=str2double(genes{k}(6:end));
end
save ~/projects/pol2rnaseq/data/rna_new_data4.mat normcounts normsds genes geneids
