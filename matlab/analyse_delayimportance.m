%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
mybasedir_code='~/synergy_data/tempcodebranch/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
mybasedir_data='~/synergy_data/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
analysisdir=[mybasedir_analyses 'analyses/'];

% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/jaakko_testversion/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/jaakko_testversion/'];
% for optimiDefaultConstraint.m
path3=[mybasedir_code 'optimi/matlab/jaakko_testversion/'];
% for lnDiffErfs.m
path4=[mybasedir_code 'ndlutil/matlab/'];
% for addPrior.m
path5=[mybasedir_code 'prior/matlab/'];
% for dist2.m
path6=[mybasedir_code 'matlab/netlab/NETLAB3p3/'];
% for modelTieParam.m
path7=[mybasedir_code 'mltools/matlab/'];
% for various experiment things
path8=[mybasedir_code 'pol2rnaseq/matlab/'];

addpath(path1,path2,path3,path4,path5,path6,path7,path8)



%---------------------------------------------------
% Load gene data
%---------------------------------------------------
cd(h3k4me3dir)
load series_for_matti_corrected_ver1.mat  % provides pol_summaryseries,rna_summaryseries

cd(analysisdir)
load models_Pol2.mat % provides interestinggenes
templls=lls_Pol2(:,1)-lls_Pol2(:,3);
[y,I]=sort(-templls);
interestinggenes=genes_Pol2;
interestinggenes=interestinggenes(I);


%---------------------------------------------------
% Load fitting results with small and long delays
%---------------------------------------------------
cd(analysisdir)
load allresults_smalldelay.mat

allresults_ensemblids_small = allresults_ensemblids; 
allresults_geneindices_small = allresults_geneindices;
allresults_jointmodels_small = allresults_jointmodels;
allresults_jointtransforminfos_small = allresults_jointtransforminfos;
allresults_loglikelihoods_small = allresults_loglikelihoods;

load allresults_longdelay.mat

allresults_ensemblids_long = allresults_ensemblids; 
allresults_geneindices_long = allresults_geneindices;
allresults_jointmodels_long = allresults_jointmodels;
allresults_jointtransforminfos_long = allresults_jointtransforminfos;
allresults_loglikelihoods_long = allresults_loglikelihoods;



%---------------------------------------------------
% Compute estimated probability of gene having long delay
%---------------------------------------------------

longdelayprobs = nan*ones(max([length(allresults_jointtransforminfos_small) ...
     length(allresults_jointtransforminfos_long)]),1);
for i=1:length(longdelayprobs),
  if (i<=size(allresults_loglikelihoods_small,1)) ...
    && (i<=size(allresults_loglikelihoods_long,1)),
    temp1=min([allresults_loglikelihoods_long(i,3) ...
	       allresults_loglikelihoods_small(i,3)]);
    if (~isempty(allresults_jointmodels_small{i})) && (~isempty(allresults_jointmodels_long{i})),
      longdelayprobs(i)=exp(allresults_loglikelihoods_long(i,3)-temp1) / ...
	  (exp(allresults_loglikelihoods_long(i,3)-temp1)+exp(allresults_loglikelihoods_small(i,3)-temp1));    
    end;
  end;  
end;


significantgenes=find(longdelayprobs>0.99);
significantgenes=interestinggenes(significantgenes);


%---------------------------------------------------
% Plot the series
%---------------------------------------------------
h1=figure;
h2=figure;
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

