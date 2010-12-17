function overlappingpeaks=find_overlapping_peaksovertime(chr_index,allpeaks,peakextentmultiplier);

overlappingpeaks=cell(1,3);

foundpeakstarts;
foundpeakends;
foundpeakheights;

n_timepoints=length(allpeaks);

peaksused=cell(n_timepoints,1);

for t=1:n_timepoints,
  peaks=allpeaks{t};
  npeaks=length(peaks{chr_index,1});

  % peaksused{t}=zeros(npeaks,1);

  for i=1:npeaks,    
    peakstart=peaks{chr_index,1}(i);
    peakend=peaks{chr_index,1}(i);
    peakcenter=(peakstart+peakend)/2;

    % extend peak according to multiplier from center
    % (note: an alternative would be to extend from summit)
    peakstart=peakcenter-peakextentmultiplier*(peakcenter-peakstart);
    peakend=peakcenter+peakextentmultiplier*(peakend-peakcenter);

    overlappingpeak=cell(n_timepoints,1);
    overlappingpeak{t}=[i];

    % find peaks in the other time points that overlap this one
    for t2=1:n_timepoints,
      if t2~=t1,
	peaks2=allpeaks{t2};
	npeaks2=length(peaks2{chr_index,index_peakstart});

	overlappingpeak{t2}=[];
	for i2=1:npeaks2,    
	  peakstart2=peaks2{chr_index,index_peakstart}(i2);
	  peakend2=peaks{chr_index,index_peakend}(i2);
	  peakcenter2=(peakstart2+peakend2)/2;

	  % extend peak according to multiplier from center
	  % (note: an alternative would be to extend from summit)
	  peakstart2=peakcenter2-peakextentmultiplier*(peakcenter2-peakstart2);
	  peakend2=peakcenter2+peakextentmultiplier*(peakend2-peakcenter2);

	  % check overlap
	  found_overlap=0;
	  if ((peakstart>=peakstart2) && (peakstart<=peakend2))
	    found_overlap=1;
	  elseif ((peakend>=peakstart2) && (peakend<=peakend2))
	    found_overlap=1;
	  elseif ((peakstart2>=peakstart) && (peakstart2<=peakend))
	    found_overlap=1;
	  end;
	  if found_overlap==1,
	    overlappingpeak{t2}=[overlappingpeak{t2} i2];
	  end;
	end;
      end;
    end;
    
    % Give different types of scores to the overlapping peak
    % -degree of overlap of peaks between time points
    %  (if there are several peaks for some time point, take the 
    %  lowest start point and highest endpoint)
    % -sum of scores (highest-scoring peak for each timepoint, if
    %  there are several)
    % -mean peak height over time
    % -variance of peak height over time
    
    % degree of overlap
    max_startpoint=-inf;
    min_endpoint=inf;
    for t2=1:n_timepoints,      
      peaks2=allpeaks{t2};
      startpoint=min(peaks{chr_index,index_peakstart}(overlappingpeaks{t2});
      endpoint=max(peaks{chr_index,index_peakend}(overlappingpeaks{t2});
      if startpoint>max_startpoint, max_startpoint=startpoint, end;
      if endpoint<min_endpoint, min_endpoint=endpoint, end;
    end;
    
    
  end;
  
end;





