function plotsinglepair(pairs,i,normalizeflag);

pol2overall_uniquehits = [101720334 91131727 72832980 93776268 93180719 94068955 87257588 85438329 92506713 94888353];
pol2filtered_hits = [82716552 76680010 51009272 78425211 74768277 73970821 59797272 56458506 83618025 81174490];
pol2filtered_multipliers = pol2filtered_hits/mean(pol2filtered_hits);

rna_uniquehits = [21695936 20435564 22212154 21769755 23547458 22137760 20776819];
rna_multipliers = rna_uniquehits/mean(rna_uniquehits);

  pol2peaks = pairs{i}{2};
  pol2peakmean = pol2peaks{10};
  pol2peakvar = pol2peaks{11};
  
  pol2peakheights = pol2peaks{9}(1:7);

  rna = cell2mat(pairs{i}(7:13));

if normalizeflag==1,
rna=rna./rna_multipliers(1:7);
end;

  clf;
  subplot(2,1,1);
  h=plot(pol2peakheights,'r-');
  set(h,'LineWidth',2);

  min_startpoint = pol2peaks{3};
  max_startpoint = pol2peaks{4};
  min_endpoint = pol2peaks{5};
  max_endpoint = pol2peaks{6};
  chrtitle=pairs{i}{4};
  
  pol2title = sprintf('POL2 at %s, overlapping peak %d - %d', chrtitle, max_startpoint, min_endpoint);
  title(pol2title);
  
  subplot(2,1,2);
  h=plot(rna,'b-');
  set(h,'LineWidth',2);
  chrtitle=pairs{i}{4};
  rnatitle=pairs{i}{3};

 % tempcorr=corr([pol2peakheights' rna']);  
 % pearsoncorr=tempcorr(1,2);
if (var(pol2peakheights)>0) && (var(rna)>0), 
pearsoncorr=((pol2peakheights-mean(pol2peakheights))*(rna-mean(rna))')/7/sqrt(var(pol2peakheights)*var(rna));
else 
pearsoncorr=0;
end;


% compute average squared deviation of RNA from zero
  rnanonzero=find(rna~=0);
  if length(rnanonzero) > 0,
    % assume that all exactly zero scores are missing values, compensate
    rnavar = sum(rna(rnanonzero).^2)*7/length(rnanonzero);
  else
    rnavar=1;
  end;
   
  % compute average squared deviation of POL2 from zero
  pol2nonzero=find(pol2peakheights~=0);
  if length(pol2nonzero) > 0,
    % assume that all exactly zero scores are missing values, compensate
    pol2var = sum(pol2peakheights(pol2nonzero).^2)*7/length(pol2nonzero);
  else
    pol2var = 1;
  end;

  % compute unnormalized correlation (just sum of products)
  sum_rna_to_pol2 = pol2peakheights*rna';

  % compute the "assumed zero-mean, missing value compensated" correlation
  zeromeancorrelation = sum_rna_to_pol2/sqrt(pol2var*rnavar);

  
  title2=sprintf('RNA at %s, %s\nPearson correlation %f, zero-mean correlation %f',...
		 chrtitle,rnatitle,pearsoncorr, zeromeancorrelation);

  title(title2);
  

