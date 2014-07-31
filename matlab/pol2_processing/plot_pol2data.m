function plot_pol2data(chr_index,mappings,peaks,bins,binsstart,interval,min_peakscore,max_peakscore,use_logscale);

if ~isempty(mappings),
  nbins=length(bins);
  bin_positions=(binsstart+interval/2)+[0:(nbins-1)]*interval;
  if use_logscale==1,
    h=plot(bin_positions,log(bins/interval+1),'b-');
  else
    h=plot(bin_positions,bins/interval,'b-');
  end;
  %set(h,'Color',[0.7 0.7 1]);
end;

hold on;


if ~isempty(peaks),
index_peakstarts=1;
index_peakends=2;
index_peakscores=3;
index_summitstarts=4;
index_summitends=5;
index_summitheights=6;

peakheight_multiplier=1;

npeaks=peaks{chr_index,8};
for i=1:npeaks,
  peakscore=peaks{chr_index,index_peakscores}(i);

  % truncate peak score to desired minimum and maximum
  if (peakscore < min_peakscore), peakscore = min_peakscore; end;
  if (peakscore > max_peakscore), peakscore = max_peakscore; end;
  propscore=(peakscore-min_peakscore)/(max_peakscore-min_peakscore);

  % use a logarithmic scale for the coloring
  peakcolor=log(propscore+1)/log(2);

  % plot peak
  h=plot([peaks{chr_index,index_peakstarts}(i) peaks{chr_index,index_peakends}(i)],[-1 -1]);
  set(h,'Color',[peakcolor 1-peakcolor 0]);

  peakheight=peaks{chr_index,index_summitheights}(i);
  peakheight=peakheight*peakheight_multiplier;

  % plot summit
  if use_logscale==1,
    h=plot([peaks{chr_index,index_summitstarts}(i) peaks{chr_index,index_summitstarts}(i)],[0 log(peakheight+1)]);
  else
    h=plot([peaks{chr_index,index_summitstarts}(i) peaks{chr_index,index_summitstarts}(i)],[0 peakheight]);
  end;
  set(h,'Color',[peakcolor 1-peakcolor 0]);
end;
end;


