% Removes duplicate pol2 peak-series from the pairs list.
% Jaakko Peltonen, Jan 11 2011.


pol2series=zeros(length(pairs),7);
for i=1:size(pairs,2),
  pol2series(i,:)=pairs{i}{2}{9}(1:7);
end;

rnaseries=zeros(length(pairs),7);
for i=1:size(pairs,2),
  rnaseries(i,:)=cell2mat(pairs{i}(7:13));
end;

uniques = ones(size(pairs,2),1);
nuniques=0;

for i=1:size(pairs,2),
  if mod(i,100)==0, i, end;
    
  duplicate_found=0;
  for j=1:nuniques,
    tempdiff=sum((pol2series(i,:)-pol2series(j,:)).^2);
    tempdiff2=sum((rnaseries(i,:)-rnaseries(j,:)).^2);
    if (tempdiff==0) && (tempdiff2==0),
      duplicate_found = 1;
      break;
    end;
  end;
  if duplicate_found==0,
    nuniques=nuniques+1;
    uniques(nuniques)=i;
  end;    
end;

uniquepairs = {pairs{uniques(1:nuniques)}};
save uniquepairs.mat uniquepairs


% Old code, more strict duplicate removal (commented out)
% Each pol2-peak-series has several peaks for each time point.
% Two pol2-peak-series are considered duplicates if they share
% any of their peaks at any of the time points.

if 1==0,

for (i=1:size(pairs,2)),
if mod(i,1000)==0, i, end;

  pol2peaks = pairs{i}{2};
  idlist=pol2peaks{1};

  duplicate_found=0;

  j=1;
  while j < length(uniques),
  
    temppeaks=pairs{uniques(j)}{2};
    idlist2=temppeaks{1};

    t=1;
    while t <= 10,
      tempids1=idlist{t};
      tempids2=idlist2{t};

      k=1;
      while k<=length(tempids1),
        k2=1;
	while k2<length(tempids2),
	  if tempids1(k)==tempids2(k2),
	    duplicate_found=1;
	    k=100000; k2=100000; t=10000; j=inf;
	  end;
	  k2=k2+1;
	end;
	k=k+1;
      end;
      t=t+1;
    end;
    j=j+1;
  end;
  
  if duplicate_found==0,
    uniques=[uniques i];
  end;
end;

end;
