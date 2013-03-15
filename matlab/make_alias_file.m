t = importdata('~/projects/pol2rnaseq/data/ensembl69_aliases_filtered.txt', '\t');
aliases = struct();
for k=1:length(t),
  aliases.(t{k}(1:15)) = t{k}(17:end);
end
