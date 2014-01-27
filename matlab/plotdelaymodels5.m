%-------------------------------------------
% Preliminaries - set up the directory paths
%-------------------------------------------

% base directory for code: set this to whatever is the correct path.
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';

% base directory for data: set this to where the result files are located.
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/analyses/';
%mybasedir_data='/share/mi/workspace/jtpelto/synergy/synergy_data/analyses/';
mybasedir_data='/share/work/jtpelto/synergy-data/analyses/';

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
load allresults_shifted_longerdelay.mat


%---------------------------------------------------
% Compute estimated probability of gene having a long-delay effect
% compared to a naive model 
%---------------------------------------------------

longdelayprobs = nan*ones(length(allresults_jointtransforminfos),1);
for i=1:length(longdelayprobs),
  if (i<=size(allresults_loglikelihoods,1)) ...
    && (i<=size(allresults_loglikelihoods,1)),
    if (~isempty(allresults_jointmodels{i})),

      % prob. of effect with long delay
      temp1=min([allresults_loglikelihoods(i,1) allresults_loglikelihoods(i,3)]);
      peffect=exp(allresults_loglikelihoods(i,3)-temp1)/...
              (exp(allresults_loglikelihoods(i,1)-temp1)+exp(allresults_loglikelihoods(i,3)-temp1));
      
      longdelayprobs(i)=peffect;
    end;
  end;  
end;

[y,I]=sort(-longdelayprobs);





timeshift=300;

%---------------------------------------------------
% Plot the series
%---------------------------------------------------
maxk=length(allresults_jointmodels);
k=1041;while k<maxk,
  k
  model=allresults_jointmodels{k};
  if ~isempty(model),
    % The next three lines are needed to fix long variable names
    % which Octave does not save correctly.
    model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
    model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
    model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

    ensemblid=allresults_ensemblids(k);
    loglik=allresults_loglikelihoods(k,3);
    plottitle=makeplottitle(model,loglik,ensemblid);
  
    predicttimes=timeshift + 1280*(([0:256]'/256).^2);
    plotpredictions(model,predicttimes,2,1,0,1,1,0,[],min(model.t));
  
    h=gca;
    ylim=get(h,'ylim');
    xlim=get(h,'xlim');
    xloc=xlim(1);
    %yloc=ylim(1)+(ylim(2)-ylim(1))*2.6238;
    yloc=ylim(1)+(ylim(2)-ylim(1))*2.4238;
    h=text(xloc,yloc,plottitle);
    set(h,'fontname','Helvetica');
    drawnow;
  
    % Octave version
    %print(sprintf('ENSG%011d_GP_shifted_delay%f.png',ensemblid,model.delay),'-depsc','-S700,600');

    % Matlab version
    h=gcf;
    set(h,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
    print(sprintf('ENSG%011d_GP_shifted_delay%f.png',ensemblid,model.delay),'-dpng','-r100'); 
  end;

  k=k+1; 
end;


if 1,
%-------------------------------------------
% write the list of genes
%-------------------------------------------
write_genelist=1;
maxk=length(allresults_jointmodels);
if write_genelist==1,
  f = fopen('largedelays_shifted_ascii.txt','w');
  for k=1:maxk,
    model=allresults_jointmodels{k};
    if ~isempty(model),
      tempdelay=allresults_jointmodels{k}.delay(1);

      temp1=min([allresults_loglikelihoods(k,1) allresults_loglikelihoods(k,3)]);
      peffect=exp(allresults_loglikelihoods(k,3)-temp1)/(exp(allresults_loglikelihoods(k,1)-temp1)+exp(allresults_loglikelihoods(k,3)-temp1));

      fprintf(f,'%d %f %e %e %f\n',...
              allresults_ensemblids(k),...
              tempdelay,...
              allresults_loglikelihoods(k,1),...
              allresults_loglikelihoods(k,3),...
              peffect); 
    end;
  end;
  fclose(f);
end;
end;
