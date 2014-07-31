pol_summaryseries_means=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
  pol_summaryseries_means(k)=mean(pol_summaryseries(k,:));
end;
pol_summaryseries_variances=zeros(size(bininfo,1),1);
for k=1:size(bininfo,1),
  pol_summaryseries_variances(k)=var(pol_summaryseries(k,:));
end;
