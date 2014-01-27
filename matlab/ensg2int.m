function id = ensg2int(gene),

if ~iscell(gene),
  id = str2double(gene(5:end));
else
  id = str2double(cellfun(@(x)(x(5:end)), gene, 'UniformOutput', false));
end
