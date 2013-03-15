function save_rna_means(r),
  
f = fopen('rna_means.txt', 'w');
for k=1:size(r.mu, 1),
  fprintf(f, '%s %f %f %f %f %f %f %f %f %f %f\n', r.geneID{k}, r.mu(k,:));
end
fclose(f);
