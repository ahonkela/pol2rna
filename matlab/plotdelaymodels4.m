%-------------------------------------------
% Preliminaries - set up the directory paths
%-------------------------------------------

% base directory for code: set this to whatever is the correct path.
mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';

% base directory for data: set this to where the result files are located.
mybasedir_data='~/jaakkos_files/synergy/synergy_data/analyses/';
%mybasedir_data='/share/mi/workspace/jtpelto/synergy/synergy_data/analyses/';

% base directory for plots: assumes the plots will be printed
% to the same directory where the data is.
mybasedir_plots=mybasedir_data;


% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/'];
% for parameter transformations etc.
path3=[mybasedir_code 'optimi/matlab/'];
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
load smallvslargedelays.mat

sortedresults_jointmodels_longdelay = sortedresults_jointmodels_longd;
sortedresults_jointmodels_smalldelay = sortedresults_jointmodels_small;
sortedresults_jointtransforminfos_longdelay = sortedresults_jointtransforminf;

%---------------------------------------------------
% Compute estimated probability of gene having a small-delay effect
% compared to a naive model 
%---------------------------------------------------


smalldelayprobs = nan*ones(max([length(sortedresults_jointtransforminfos_longdelay) ...
     length(sortedresults_jointtransforminfos_longdelay)]),1);
for i=1:length(smalldelayprobs),
  if (i<=size(sortedresults_loglikelihoods,1)) ...
    && (i<=size(sortedresults_loglikelihoods,1)),
    if (~isempty(sortedresults_jointmodels_smalldelay{i})) && (~isempty(sortedresults_jointmodels_longdelay{i})),

      % prob. of effect with small delay
      temp1=min([sortedresults_loglikelihoods(i,1) sortedresults_loglikelihoods(i,2)]);
      peffect=exp(sortedresults_loglikelihoods(i,2)-temp1)/...
              (exp(sortedresults_loglikelihoods(i,1)-temp1)+exp(sortedresults_loglikelihoods(i,2)-temp1));
      
      smalldelayprobs(i)=peffect;    
    end;
  end;  
end;

[y,I]=sort(-smalldelayprobs);

k=1;
  model=sortedresults_jointmodels_smalldelay{I(k)};
  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
plottitle='';
predicttimes=1280*(([0:256]'/256).^2);  
plotpredictions(model,predicttimes,2,1,1,plottitle);
drawnow;

eid=??;
k=find(sortedresults_ensemblids==eid);
  model=sortedresults_jointmodels_longdelay{k};
  %model=sortedresults_jointmodels_smalldelay{k};
  model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;
  ensemblid=sortedresults_ensemblids(k);
  loglik=sortedresults_loglikelihoods(k,3);
  %loglik=sortedresults_loglikelihoods(k,2);
  plottitle=makeplottitle(model,loglik,ensemblid);
predicttimes=1280*(([0:256]'/256).^2);  
plotpredictions(model,predicttimes,2,1,1,[]);
h=gca;
ylim=get(h,'ylim');
yloc=ylim(1)+(ylim(2)-ylim(1))*4.0857;
text(0,yloc,plottitle);
drawnow;
print(sprintf('ENSG%011d_GP_delay%f.eps',sortedresults_ensemblids(k),model.delay),'-depsc','-S800,1000');



%---------------------------------------------------
% Plot the series
%---------------------------------------------------
maxk=length(sortedresults_jointmodels_longd);
k=1;while k<maxk,
  k
  model=sortedresults_jointmodels_longd{k};
  % The next three lines are needed to fix long variable names
  % which Octave does not save correctly.
  model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
  model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
  model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

  % title of the plot
  temp1=min([sortedresults_loglikelihoods(k,1) sortedresults_loglikelihoods(k,3)]);
  peffect=exp(sortedresults_loglikelihoods(k,3)-temp1)/(exp(sortedresults_loglikelihoods(k,1)-temp1)+exp(sortedresults_loglikelihoods(k,3)-temp1));
  temp1=min([sortedresults_loglikelihoods(k,2) sortedresults_loglikelihoods(k,3)]);
  plongdelay=exp(sortedresults_loglikelihoods(k,3)-temp1)/(exp(sortedresults_loglikelihoods(k,2)-temp1)+exp(sortedresults_loglikelihoods(k,3)-temp1));
  
  plottitle=...
      sprintf('Gene %d (ENSEMBL %d), loglik %f, pEffect=%f, pLongDelay=%f\nB=%f, D=%f, S=%f, delay=%f, MuPol=%f, Rna0=%f, nvarPol=%f\n',...
	      k, ...                                % running index of the gene
	      sortedresults_ensemblids(k), ...            % ensembl-id of the gene
	      sortedresults_loglikelihoods(k,3), ...      % log-likelihood of the long-delay model
              peffect,plongdelay,...
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
  drawnow;
      
  % save the plot as an EPS file
  printplots=1;
  if printplots==1,
    % ENSG00000196208_GP.png
    print(sprintf('ENSG%011d_GP.png',sortedresults_ensemblids(k)),'-dpng','-S400,400');
  end;
  close all
 
  k=k+1; 
%  x=kbhit();
%  if x=='w',k=k-1;
%  else k=k+1;end;
end;

if 0,
%-------------------------------------------
% write a LaTeX file that will display the plots
%-------------------------------------------
write_latexfile=1;
if write_latexfile==1,
  f = fopen('delayplots.tex','w');
  fprintf(f,'\\documentclass{article}\n')
  fprintf(f,'\\usepackage{graphicx}\n')
  fprintf(f,'\\addtolength{\\textwidth}{11cm}\n');
  fprintf(f,'\\addtolength{\\hoffset}{-4.5cm}\n');
  fprintf(f,'\\begin{document}\n');
  for k=1:maxk,
    if k>1, fprintf(f,'\\newpage\n'); end;
    fprintf(f,'$\\phantom{a}$\n');
    fprintf(f,'\\hspace{-2cm}\n');
    fprintf(f,'\\includegraphics[width=\\textwidth]{tempprint%d.eps}\n',k);  
  end;
  fprintf(f,'\\end{document}');
  fclose(f);
end;
end;

if 0,
%-------------------------------------------
% write the list of genes
%-------------------------------------------
write_genelist=1;
maxk=length(sortedresults_jointmodels_longd);
if write_genelist==1,
  f = fopen('smallvslargedelays_ascii.txt','w');
  for k=1:maxk,
    tempdelay=sortedresults_jointmodels_longd{k}.delay(1);

    temp1=min([sortedresults_loglikelihoods(k,1) sortedresults_loglikelihoods(k,3)]);
    peffect=exp(sortedresults_loglikelihoods(k,3)-temp1)/(exp(sortedresults_loglikelihoods(k,1)-temp1)+exp(sortedresults_loglikelihoods(k,3)-temp1));

    temp1=min([sortedresults_loglikelihoods(k,2) sortedresults_loglikelihoods(k,3)]);
    plongdelay=exp(sortedresults_loglikelihoods(k,3)-temp1)/(exp(sortedresults_loglikelihoods(k,2)-temp1)+exp(sortedresults_loglikelihoods(k,3)-temp1));

    fprintf(f,'%d %f %e %e %e %f %f\n',...
            sortedresults_ensemblids(k),...
            tempdelay,...
            sortedresults_loglikelihoods(k,1),...
            sortedresults_loglikelihoods(k,2),...
            sortedresults_loglikelihoods(k,3),...
            peffect, plongdelay); 
  end;
  fclose(f);
end;
end;
