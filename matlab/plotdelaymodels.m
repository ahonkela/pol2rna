%---------------------------------------------------
% Plot the series
%---------------------------------------------------
k=1;while k<length(significantgenes),
  gene_index=significantgenes(k);
  tempseries=pol_summaryseries(gene_index,:);
  tempseries=(tempseries-min(tempseries))./(max(tempseries)-min(tempseries));

  model=allresults_jointmodels_small{gene_index};  
  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

  % title of the plot
  plottitle=...
      sprintf('Gene %d (ENSEMBL %d), loglik %f. B=%f, D=%f, S=%f,\ndelay=%f, MuPol=%f, Rna0=%f, nvarPol=%f\n',...
	      k, ...                                % running index of the gene
	      allresults_ensemblids_small(k), ...            % ensembl-id of the gene
	      allresults_loglikelihoods_small(k,3), ...      % log-likelihood of the model
	      model.B(1), ...                       % basal rate of RNA production
	      model.D(1), ...                       % RNA decay rate
	      model.S(1), ...                       % RNA sensitivity to PolII
	      model.delay(1), ...                   % delay between PolII and RNA
	      model.simMean, ...                    % PolII mean (often closest to first observations)
	      model.disimStartMean, ...             % initial RNA concentration
	      model.kern.comp{2}.comp{1}.variance); % Observation-noise variance for PolII     
    
  predicttimes=1280*(([0:256]'/256).^2);  
  % create the plot
  figure;
  plotpredictions(model,predicttimes,2,1,1,plottitle);
  

  model=allresults_jointmodels_long{gene_index};  
  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

  % title of the plot
  plottitle=...
      sprintf('Gene %d (ENSEMBL %d), loglik %f. B=%f, D=%f, S=%f,\ndelay=%f, MuPol=%f, Rna0=%f, nvarPol=%f\n',...
	      k, ...                                % running index of the gene
	      allresults_ensemblids_long(k), ...            % ensembl-id of the gene
	      allresults_loglikelihoods_long(k,3), ...      % log-likelihood of the model
	      model.B(1), ...                       % basal rate of RNA production
	      model.D(1), ...                       % RNA decay rate
	      model.S(1), ...                       % RNA sensitivity to PolII
	      model.delay(1), ...                   % delay between PolII and RNA
	      model.simMean, ...                    % PolII mean (often closest to first observations)
	      model.disimStartMean, ...             % initial RNA concentration
	      model.kern.comp{2}.comp{1}.variance); % Observation-noise variance for PolII     
    
  predicttimes=1280*(([0:256]'/256).^2);  
  % create the plot
  figure;
  plotpredictions(model,predicttimes,2,1,1,plottitle);
  
  
%  plot(sqrt([0 5 10 20 40 80 160 320 640 1280]), ...
%       tempseries,'r-');
%  title(sprintf('%d',k));
%  hold on;
%  tempseries=rna_summaryseries(gene_index,:);
%  tempseries=(tempseries-min(tempseries))./(max(tempseries)-min(tempseries));
%  plot(sqrt([0 5 10 20 40 80 160 320 640 1280]), ...
%       tempseries,'g-');  
  
  x=kbhit();
  if x=='w',k=k-1;
  else k=k+1;end;
end;

