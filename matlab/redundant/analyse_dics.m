r = {};
for k=1:16, r{k} = load(sprintf('dics_2012-01-23_%d.mat', k)); end
dic_joint = r{1}.dic_joint; for k=2:16, dic_joint = dic_joint + r{k}.dic_joint; end
dic_pol = r{1}.dic_pol; for k=2:16, dic_pol = dic_pol + r{k}.dic_pol; end
dic_rna = r{1}.dic_rna; for k=2:16, dic_rna = dic_rna + r{k}.dic_rna; end
