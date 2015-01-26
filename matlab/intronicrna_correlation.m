
load ~/projects/pol2rnaseq/data/rpkm.mat

load ~/projects/pol2rnaseq/data/bininfo_nov2014_corrected.mat
load ~/projects/pol2rnaseq/data/pol2_summaryseries_2014_11_19.mat

% Ensembl ids for Hande's genes as numbers
hande_geneids = nan*ones(length(genes),1);
for k=1:length(genes),
  hande_geneids(k) = str2double(genes{k}(5:end));
end;

pol2intronicrna_correlation = nan*ones(size(bininfo,1),1);
for k=1:size(bininfo,1),
  % Find the gene in Hande's array
  k2 = find(hande_geneids==bininfo(k,5));
  if length(k2==1),
    pol2series=pol2_summaryseries(k,:);
    intronicrnaseries=mT2(k2,:);
    tempcorr = corr(pol2series',intronicrnaseries');
    pol2intronicrna_correlation(k)=tempcorr;
  end;
end;

save pol2intronicrna_correlation.mat bininfo pol2intronicrna_correlation -mat

f=fopen('pol2intronicrna_correlation.txt','w');
for k=1:size(bininfo,1),
  if isnan(pol2intronicrna_correlation(k))==0,
    fprintf(f,'ENSG%011d %f\n', bininfo(k,5), pol2intronicrna_correlation(k));
  end;
end;
fclose(f);

