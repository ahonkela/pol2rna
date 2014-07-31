%mybasedir_code='/share/work/jtpelto/tempsynergy/';
%mybasedir_code='/media/JPELTONEN4/mlprojects/';
mybasedir_code='/share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/';
%mybasedir_code='~/jaakkos_files/synergy/mlprojects/';
%mybasedir_code='~/mlprojects/';
%mybasedir_code='/share/work/jtpelto/tempsynergy/';

%mybasedir_data='/share/work/jtpelto/tempsynergy/';
%mybasedir_data='/media/JPELTONEN4/synergy_data/';
%mybasedir_data='~/jaakkos_files/synergy/synergy_data/';
%mybasedir_data='~/synergy_data/';
mybasedir_data='/share/work/jtpelto/synergy-data/';

mybasedir_analyses=mybasedir_data;
%mybasedir_analyses='~/jaakkos_files/synergy/';


% pol2 and H3K4me3 data
pol2dir=[mybasedir_data 'PolII/Mapping_results/'];
h3k4me3dir=[mybasedir_data 'H3K4me3/Mapping_results/'];
analysisdir=[mybasedir_analyses 'analyses/'];

% for kernel-level computations
path1=[mybasedir_code 'kern/matlab/']
% for model-level computations
path2=[mybasedir_code 'gpsim/matlab/'];
% for optimiDefaultConstraint.m
path3=[mybasedir_code 'optimi/matlab/'];
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
% cd(h3k4me3dir)
% load series_for_matti_corrected_ver1.mat  % provides pol_summaryseries,rna_summaryseries

% cd(analysisdir)
% load models_Pol2.mat % provides interestinggenes
% templls=lls_Pol2(:,1)-lls_Pol2(:,3);
% [y,I]=sort(-templls);
% interestinggenes=genes_Pol2;
% interestinggenes=interestinggenes(I);


%---------------------------------------------------
% Load fitting results with small and long delays
%---------------------------------------------------
cd(analysisdir)

load allresults_shifted_polorrna4.mat
% allresults_ensemblids_pol2 allresults_geneindices_pol2 allresults_jointmodels_pol2 allresults_jointtransforminfos_pol2 allresults_loglikelihoods_pol2
% allresults_ensemblids_rna allresults_geneindices_rna allresults_jointmodels_rna allresults_jointtransforminfos_rna allresults_loglikelihoods_rna
allresults_jointtransforminfos_pol2 = allresults_tinfos_pol2;
allresults_jointtransforminfos_rna = allresults_tinfos_rna;


load allresults_shifted_longerdelay5.mat
allresults_ensemblids_joint = allresults_ensemblids; 
allresults_geneindices_joint = allresults_geneindices;
allresults_jointmodels_joint = allresults_jointmodels;
allresults_jointtransforminfos_joint = allresults_jointtransforminfos;
allresults_loglikelihoods_joint = allresults_loglikelihoods;




%---------------------------------------------------
% Compute estimated probability of gene having long delay
%---------------------------------------------------

probcomparison = nan*ones(max([length(allresults_jointtransforminfos_pol2) ...
     length(allresults_jointtransforminfos_joint)]),1);
for i=1:length(probcomparison),
  if (i<=size(allresults_loglikelihoods_joint,1)) ...
    && (i<=size(allresults_loglikelihoods_pol2,1)),
    if (~isempty(allresults_jointmodels_joint{i})) && (~isempty(allresults_jointmodels_pol2{i})),
      probcomparison(i)=allresults_loglikelihoods_joint(i,3) ...
          -allresults_loglikelihoods_pol2(i,3) -allresults_loglikelihoods_rna(i,3);
    end;
  end;  
end;


[y,I]=sort(-probcomparison);







%---------------------------------------------------
% Plot the series
%---------------------------------------------------
timeshift=300;

maxl=sum(isnan(y)==0);
%maxl=100;
%l=1;
l=5836;
while l<maxl,
  k=I(l)
  model=allresults_jointmodels_joint{k};
  if ~isempty(model),
    % The next three lines are needed to fix long variable names
    % which Octave does not save correctly.
    model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
    model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
    model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

    ensemblid=allresults_ensemblids_joint(k);
    loglik=allresults_loglikelihoods_joint(k,3);
    plottitle=makeplottitle(model,loglik,ensemblid);
    plottitle=['JOINTFIT ' plottitle];

    predicttimes=timeshift + 1280*(([0:256]'/256).^2);
    plotpredictions(model,predicttimes,2,1,0,1,1,0,[],timeshift);
  
%    h=gca;
%    ylim=get(h,'ylim');
%    xlim=get(h,'xlim');
%    xloc=xlim(1);
%    %yloc=ylim(1)+(ylim(2)-ylim(1))*2.6238;
%    %yloc=ylim(1)+(ylim(2)-ylim(1))*2.4238;
%    yloc=ylim(1)+(ylim(2)-ylim(1))*8.4238;
%    h=text(xloc,yloc,plottitle);
%    set(h,'fontname','Helvetica');
%    h=gca;
%    set(h,'ylim',ylim);
%    set(h,'xlim',xlim);
%    drawnow;
    axes;
    axis([0 1 0 1]);
    h=text(0,1.07,plottitle);
    %h=title(plottitle);
    set(h,'fontname','Helvetica');    
    set(h,'fontsize',9);
    axis off;
 
    % Octave version
    % print(sprintf('ENSG%011d_GP_shifted_delay%f.eps',ensemblid,model.delay),'-depsc','-S700,600');
    print(sprintf('ENSG%011d_gppol2rna.png',ensemblid),'-dpng','-S700,600');

    % % Matlab version
    % h=gcf;
    % set(h,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
    % %print(sprintf('ENSG%011d_GP_shifted_delay%f.png',ensemblid,model.delay),'-dpng','-r100'); 
  end;


  model=allresults_jointmodels_pol2{k};
  if ~isempty(model),
    % The next three lines are needed to fix long variable names
    % which Octave does not save correctly.
    model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
    model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
    model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

    ensemblid=allresults_ensemblids_pol2(k);
    loglik=allresults_loglikelihoods_pol2(k,3);
    plottitle=makeplottitle(model,loglik,ensemblid);
    plottitle=['POL2FIT ' plottitle];
  
    predicttimes=timeshift + 1280*(([0:256]'/256).^2);
    plotpredictions(model,predicttimes,2,1,0,1,1,0,[],timeshift);
  
    axes;
    axis([0 1 0 1]);
    h=text(0,1.07,plottitle);
    set(h,'fontname','Helvetica');    
    set(h,'fontsize',9);
    axis off;
 
    % Octave version
    % print(sprintf('ENSG%011d_GP_shifted_delay%f_pol2.eps',ensemblid,model.delay),'-depsc','-S700,600');
    print(sprintf('ENSG%011d_gppol2.png',ensemblid),'-dpng','-S700,600');
  end;

  model=allresults_jointmodels_rna{k};
  if ~isempty(model),
    % The next three lines are needed to fix long variable names
    % which Octave does not save correctly.
    model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
    model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
    model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

    ensemblid=allresults_ensemblids_rna(k);
    loglik=allresults_loglikelihoods_rna(k,3);
    plottitle=makeplottitle(model,loglik,ensemblid);
    plottitle=['RNAFIT ' plottitle];
  
    predicttimes=timeshift + 1280*(([0:256]'/256).^2);
    plotpredictions(model,predicttimes,2,1,0,1,1,0,[],timeshift);
  
    axes;
    axis([0 1 0 1]);
    h=text(0,1.07,plottitle);
    set(h,'fontname','Helvetica');    
    set(h,'fontsize',9);
    axis off;
 
    % Octave version
    % print(sprintf('ENSG%011d_GP_shifted_delay%f_rna.eps',ensemblid,model.delay),'-depsc','-S700,600');
    print(sprintf('ENSG%011d_gprna.png',ensemblid),'-dpng','-S700,600');
  end;

  l=l+1; 
end;



if 1,
f = fopen('polrnaconnectionresults.txt','w');
fprintf(f,'ENSID LL_naive LL_indep LL_joint LL_pol2 LL_rna delay\n');
timeshift=300;
maxl=sum(isnan(y)==0);
%maxl=100;
l=1;
%l=5836;
while l<maxl,
  k=I(l)
  model=allresults_jointmodels_joint{k};
  %if ~isempty(model),
    % The next three lines are needed to fix long variable names
    % which Octave does not save correctly.
    model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
    model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
    model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

    ensemblid=allresults_ensemblids_joint(k);
    loglikjoint=allresults_loglikelihoods_joint(k,3);
    loglikpol2=allresults_loglikelihoods_pol2(k,3);
    loglikrna=allresults_loglikelihoods_rna(k,3);
    loglikindep=loglikpol2+loglikrna;
    logliknaive=allresults_loglikelihoods_joint(k,1);

  fprintf(f,'ENSG%011d %e %e %e %e %e %e\n', ensemblid, logliknaive, loglikindep, loglikjoint, loglikpol2, loglikrna, model.delay);


  %end;
  l=l+1;
end;
fclose(f);

end;


if 0,
%-------------------------------------------
% write a LaTeX file that will display the plots
%-------------------------------------------
write_latexfile=1;
if write_latexfile==1,
  f = fopen('polrnaconnectionplots.tex','w');
  fprintf(f,'\\documentclass{article}\n')
  fprintf(f,'\\usepackage{graphicx}\n')
  %fprintf(f,'\\addtolength{\\textwidth}{11cm}\n');
  fprintf(f,'\\addtolength{\\textheight}{11cm}\n');
  fprintf(f,'\\addtolength{\\voffset}{-4.5cm}\n');
  fprintf(f,'\\begin{document}\n');
  fprintf(f,'\\pagestyle{empty}\n');
  fprintf(f,'\\thispagestyle{empty}\n');
  maxl=100;
  l=1;while l<maxl,
    k=I(l)
    ensemblid=allresults_ensemblids_joint(k);
    model=allresults_jointmodels_joint{k};
    name1=sprintf('ENSG%011d_GP_shifted_delay%f.eps',ensemblid,model.delay);
    model=allresults_jointmodels_pol2{k};
    name2=sprintf('ENSG%011d_GP_shifted_delay%f_pol2.eps',ensemblid,model.delay);
    model=allresults_jointmodels_rna{k};
    name3=sprintf('ENSG%011d_GP_shifted_delay%f_rna.eps',ensemblid,model.delay);

    if l>1, fprintf(f,'\\newpage\n'); end;
    fprintf(f,'$\\phantom{a}$\\\\\n');
    fprintf(f,'\\vspace{-0.5cm}\n');
    fprintf(f,'\\includegraphics[width=0.9\\textwidth]{%s}\\\\\n',name1);
    fprintf(f,'\\vspace{-0.5cm}\n');
    fprintf(f,'\\includegraphics[width=0.9\\textwidth]{%s}\\\\\n',name2);
    fprintf(f,'\\includegraphics[width=0.9\\textwidth]{%s}\\\\\n',name3);
    l=l+1;
  end;
  fprintf(f,'\\end{document}');
  fclose(f);

end;
end;



