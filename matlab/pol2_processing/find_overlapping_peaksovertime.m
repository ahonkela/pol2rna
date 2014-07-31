function overlappingpeaks=find_overlapping_peaksovertime(chr_index,allpeaks,peakextentmultiplier);


n_timepoints=length(allpeaks);

index_peakstart=1;
index_peakend=2;
index_peakscore=3;
index_summitstart=4;
index_summitend=5;
index_summitheight=6;
index_npeaks=7;
index_maxpeaks=8;
index_linenumber=9;


overlappingpeaks={};


for t=1:n_timepoints,
  peaks=allpeaks{t};
  npeaks=length(peaks{chr_index,1});

  for i=1:npeaks,
    if mod(i,100)==0,
      fprintf(1,'Processing overlaps starting from time point %d, peak %d of %d\n',t, i, npeaks);
    end;
 
    peakstart=peaks{chr_index,index_peakstart}(i);
    peakend=peaks{chr_index,index_peakend}(i);
    peakcenter=(peakstart+peakend)/2;

    % extend peak according to multiplier from center
    % (note: an alternative would be to extend from summit)
    peakstart=peakcenter-peakextentmultiplier*(peakcenter-peakstart);
    peakend=peakcenter+peakextentmultiplier*(peakend-peakcenter);

    overlappingpeak=cell(n_timepoints,1);
    overlappingpeak{t}=[i];

    % find peaks in the other time points that overlap this one
    for t2=1:n_timepoints,
      if t2~=t,
	peaks2=allpeaks{t2};
	npeaks2=length(peaks2{chr_index,index_peakstart});

	overlappingpeak{t2}=[];
	for i2=1:npeaks2,    
	  peakstart2=peaks2{chr_index,index_peakstart}(i2);
	  peakend2=peaks2{chr_index,index_peakend}(i2);
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

    
    % extent of all peaks (union) and degree of overlap (intersection). Ignore time points without peaks.
    min_startpoint=inf;
    max_startpoint=-inf;
    min_endpoint=inf;
    max_endpoint=-inf;
    for t2=1:n_timepoints,      
      if length(overlappingpeak{t2})>0,
        peaks2=allpeaks{t2};    
	startpoint=min(peaks2{chr_index,index_peakstart}(overlappingpeak{t2}));
	endpoint=max(peaks2{chr_index,index_peakend}(overlappingpeak{t2}));
	if startpoint>max_startpoint, max_startpoint=startpoint; end;
	if startpoint<min_startpoint, min_startpoint=startpoint; end;
	if endpoint<min_endpoint, min_endpoint=endpoint; end;
	if endpoint>max_endpoint, max_endpoint=endpoint; end;
      end;      
    end;
    degree_of_overlap = min_endpoint - max_startpoint;
    peaks_extent = max_endpoint - min_startpoint;
    
    % sum of scores (time points without peaks get zero score)
    scoresum=0;
    for t2=1:n_timepoints,      
      if length(overlappingpeak{t2})>0,
	peaks2=allpeaks{t2};
	score=max(peaks2{chr_index,index_peakscore}(overlappingpeak{t2}));
	scoresum=scoresum+score;
      end;
    end;

    % mean and variance of the summit height over time; if there are several
    % peaks for some time point, use the highest summit.
    % If some time point does not have a peak, it counts as zero (this could be
    % done better using the original read data mappings; also, then we could use
    % something else than just the summit height).
    summitheights=[];
    for t2=1:n_timepoints,      
      if length(overlappingpeak{t2})>0,
	peaks2=allpeaks{t2};	
	tempheight=max(peaks2{chr_index,index_summitheight}(overlappingpeak{t2}));
	summitheights(t2)=tempheight;
      else
	summitheights(t2)=0;
      end;
    end;
    summitheight_mean=mean(summitheights);
    summitheight_variance=var(summitheights);
    

    % compute the other scores too, TODO


    % store the overlapping peak and its scores
    peakandscores = {overlappingpeak, peaks_extent, min_startpoint, max_startpoint, min_endpoint, max_endpoint, degree_of_overlap, scoresum, summitheights, summitheight_mean, summitheight_variance};
    overlappingpeaks = {overlappingpeaks{:}, peakandscores};  
    % it's inefficient to expand the array all the time like above,
    % it would be better to expand it with ever larger increments like in read_peak_summit_files.m

  end;
  
end;





