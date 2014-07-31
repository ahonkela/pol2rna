pol2overall_uniquehits = [101720334 91131727 72832980 93776268 93180719 94068955 87257588 85438329 92506713 94888353];
pol2filtered_hits = [82716552 76680010 51009272 78425211 74768277 73970821 59797272 56458506 83618025 81174490];
pol2filtered_multipliers = pol2filtered_hits/mean(pol2filtered_hits);

% Computes correlations for the time-series pairs
myscores=zeros(length(pairs),6);
for (i=1:size(pairs,2)),

  pol2peaks = pairs{i}{2};
  pol2peakmean = pol2peaks{10};
  pol2peakvar = pol2peaks{11};
  
  pol2peakheights = pol2peaks{9}(1:7);
  % pol2peakheights = pol2peakheights./pol2filtered_multipliers(1:7); % adjust for overall hit level per time point

  rna = cell2mat(pairs{i}(7:13));

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
  
  %tempcorr=corr([pol2peakheights' rna']);
  %myscores(i,1) = tempcorr(1,2);
  if (var(pol2peakheights)>0) && (var(rna)>0), 
    myscores(i,1)=(pol2peakheights-mean(pol2peakheights))*(rna-mean(rna))'/sqrt(var(pol2peakheights)*var(rna));
  else 
    myscores(i,1)=0;
  end;
  myscores(i,2) = pol2peakheights*rna';
  myscores(i,3) = pol2peakmean;
  myscores(i,4) = pol2peakvar;
  myscores(i,5) = zeromeancorrelation;

  if (length(pol2nonzero)==7) && (length(rnanonzero)==7),
    myscores(i,6)=1;
  else
    myscores(i,6)=0;
  end;
  
  if myscores(i,5)>1,
    fprintf(1, 'i=%d, myscores(i,2)=%d, pol2var=%d, rnavar=%d\n', i, myscores(i,2), pol2var,rnavar);
    pause
  end;
  
  
  % co = corr([pol2; rna]');
  %pairs{i}{1} = co(1,2); 
  if (mod(i,1000) == 0),
    fprintf('Peak %d\n', i);
  end
end

[y,I]=sort(-myscores(:,5)-abs(myscores(:,1)));

