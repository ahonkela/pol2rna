pol2overall_uniquehits = [101720334 91131727 72832980 93776268 93180719 94068955 87257588 85438329 92506713 94888353];
pol2filtered_hits = [82716552 76680010 51009272 78425211 74768277 73970821 59797272 56458506 83618025 81174490];
pol2filtered_multipliers = pol2filtered_hits/mean(pol2filtered_hits);



pol2series=zeros(length(uniquepairs),10);
for i=1:size(uniquepairs,2),
  pol2series(i,:)=uniquepairs{i}{2}{9}(1:10);
  % adjust for overall hit level
  pol2series(i,:)=pol2series(i,:)./pol2filtered_multipliers;

  % normalize
  pol2series(i,:)=pol2series(i,:)/sum(pol2series(i,:));
end;



uniqueseries=pol2series;
idx=kmeans(uniqueseries,20,'EmptyAction','singleton','MaxIter',1000);


nmembers=zeros(20,1);
for i=1:20,
I=find(idx==i);
nmembers(i)=length(I);
end;
[y,I2]=sort(-nmembers);


for k=1:5,
figure;
for j0=1:4,
subplot(2,2,j0); j=(k-1)*4+j0;

I=find(idx==I2(j));

clustmean=mean(uniqueseries(I,:));
clustvar=var(uniqueseries(I,:));

boxplot(uniqueseries(I,:));
hold on; h=plot(clustmean,'k-');
set(h,'LineWidth',2);
mytitle=sprintf('POL2 peak-column cluster %d (%d members)', j, length(I));
title(mytitle);
end;
end;




nonzeroseries=find(sum(pol2series==0,2)==0);
uniqueseries=pol2series(nonzeroseries,:);
idx=kmeans(uniqueseries,20,'EmptyAction','singleton','MaxIter',1000);


nmembers=zeros(20,1);
for i=1:20,
  I=find(idx==i);
  nmembers(i)=length(I);
end;
[y,I2]=sort(-nmembers);


for k=1:5,
  figure;
  for j0=1:4,
    subplot(2,2,j0); j=(k-1)*4+j0;

    I=find(idx==I2(j));

    clustmean=mean(uniqueseries(I,:));
    clustvar=var(uniqueseries(I,:));
    
    boxplot(uniqueseries(I,:));
    hold on; h=plot(clustmean,'k-');
    set(h,'LineWidth',2);
    mytitle=sprintf('POL2 nonzero-peak-column cluster %d (%d members)', j, length(I));
    title(mytitle);
  end;
end;




