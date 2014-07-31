I=find(isnan(sum(pol_summaryseries,2))==0);
polmean=mean(mean(pol_summaryseries(I,:)));

I=find(isnan(sum(groseq_summaryseries,2))==0);
groseqmean=mean(mean(groseq_summaryseries(I,:)));

figure;
for k=1:length(interestinggenes),
  temppol=pol_summaryseries(interestinggenes(k),:);
  tempgroseq=groseq_summaryseries(interestinggenes(k),:);
  temppol=temppol/polmean;
  tempgroseq=tempgroseq/groseqmean;
  clf;
  h=plot(sqrt([0 5 10 20 40 80 160]),temppol(1:7),'r*-');
  set(h,'LineWidth',2);
  set(h,'MarkerSize',20);
  hold on;
  h=plot(sqrt([0 0 10 10 40 40 160 160]),tempgroseq,'k*');
  set(h,'MarkerSize',20);
  pause
end;
