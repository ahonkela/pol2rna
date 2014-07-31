%-------------------------------------------
% Preliminaries - set up the directory paths
%-------------------------------------------

% base directory for code: this assumes that the current 
% matlab file is located at the root of the code directory.
tempstring=which('plot_gpnddisim_results');
mybasedir_code=tempstring(1:end-length('plot_gpnddisim_results.m'));

% base directory for data: assumes the data is stored at 
% the root of the code directory (not wise generally, this 
% is just for simple distribution of data and code together).

%mybasedir_data=mybasedir_code;
mybasedir_data='/home/jaakkopeltonen/synergy_data/analyses/';

% base directory for plots: assumes the plots will be printed
% to the same directory where the data is.
mybasedir_plots=mybasedir_data;


% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/jaakko_testversion/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/jaakko_testversion/'];
% for parameter transformations etc.
path3=[mybasedir_code 'optimi/matlab/jaakko_testversion/'];
% for lnDiffErfs.m etc.
path4=[mybasedir_code 'ndlutil/matlab/'];
% for addPrior.m etc.
path5=[mybasedir_code 'prior/matlab/'];
% for dist2.m
path6=[mybasedir_code 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m etc.
path7=[mybasedir_code 'mltools/matlab/'];
% for various experiment things and plotting code
path8=[mybasedir_code 'pol2rnaseq/matlab/'];
addpath(path1,path2,path3,path4,path5,path6,path7,path8);


%-------------------------------------------
% Load the fitted models (which include their data)
%-------------------------------------------

cd(mybasedir_data);
%load fittingresults_corrected_temp_0.mat
load allfittingresults_pol2end_rna.mat


%-------------------------------------------
% Plot each model into a file
%-------------------------------------------
drawplots=1;
printplots=0;

if drawplots==1,
  cd(mybasedir_plots);
  for k=1:length(allresults_jointmodels),
    model=allresults_jointmodels{k};
    
    if ~isempty(model),
      
      % The next three lines are needed to fix long variable names
      % which Octave does not save correctly.
      model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
      model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
      model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
      
      % title of the plot
      plottitle=...
	  sprintf('Gene %d (ENSEMBL %d), loglik %f. B=%f, D=%f, S=%f,\ndelay=%f, MuPol=%f, Rna0=%f, nvarPol=%f\n',...
		  k, ...                                % running index of the gene
		  allresults_ensemblids(k), ...            % ensembl-id of the gene
		  allresults_loglikelihoods(k,3), ...      % log-likelihood of the model
		  model.B(1), ...                       % basal rate of RNA production
		  model.D(1), ...                       % RNA decay rate
		  model.S(1), ...                       % RNA sensitivity to PolII
		  model.delay(1), ...                   % delay between PolII and RNA
		  model.simMean, ...                    % PolII mean (often closest to first observations)
		  model.disimStartMean, ...             % initial RNA concentration
		  model.kern.comp{2}.comp{1}.variance); % Observation-noise variance for PolII 
      
      % compute the posterior curves at points with time-squared
      % spacing, up to 1280 minutes
      predicttimes=1280*(([0:256]'/256).^2);
      
      % create the plot
      plotpredictions(model,predicttimes,2,1,1,plottitle);
      drawnow;
      
      % save the plot as an EPS file
      if printplots==1,
	print(sprintf('tempprint%d.eps',k),'-depsc');
      end;
      
    end;
  end;
end;


%-------------------------------------------
% write a LaTeX file that will display the plots
%-------------------------------------------
write_latexfile=0;
if write_latexfile==1,
  f = fopen('allplots.tex','w');
  fprintf(f,'\\documentclass{article}\n')
  fprintf(f,'\\usepackage{graphicx}\n')
  fprintf(f,'\\addtolength{\\textwidth}{11cm}\n');
  fprintf(f,'\\addtolength{\\hoffset}{-4.5cm}\n');
  fprintf(f,'\\begin{document}\n');
  for k=1:length(allresults_jointmodels),
    if k>1, fprintf(f,'\\newpage\n'); end;
    fprintf(f,'$\\phantom{a}$\n');
    fprintf(f,'\\hspace{-2cm}\n');
    fprintf(f,'\\includegraphics[width=\\textwidth]{tempprint%d.eps}\n',k);  
  end;
  fprintf(f,'\\end{document}');
  fclose(f);
end;


