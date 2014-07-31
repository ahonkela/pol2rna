function [d, genes] = load_bitseq_data(),
  
filepath = '~/projects/pol2rnaseq/data/gene_summary_files/';
files = {'t0000', 't0005', 't0010', 't0020', 't0040',...
         't0080', 't0160', 't0320', 't0640', 't1280'};

dc = cell(10, 1);
for k=1:length(files),
  dc{k} = importdata([filepath, files{k}]);
end
d = cat(3, dc{:});
genes = importdata([filepath, 'gene_names']);
