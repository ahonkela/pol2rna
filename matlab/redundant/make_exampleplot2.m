
allids={};
allnames={};
%Genes with reasonably early Pol2 peaks:
ids1 = [186628 249822 238117    163659    184254     170500];
names1={'FSD2','',    '',       'TIPARP', 'ALDH1A3', 'LONRF2'};
allids{1} = ids1;
allnames{1} = names1;

% Genes with fairly long delay:
ids2 = [125735    139514    142207  107554   135763 111641 133639 134897 68784   52749   131979];
names2={'TNFSF14','SLC7A1', 'URB1', 'DNMBP','URB2', 'NOP2','BTG1','BIVM','SRBD1','RRP12','GCH1'};
allids{2} = ids2;
allnames{2} = names2;

% Genes with small to medium delay:
ids3 = [  217791 230807 60303  113369   73282  152413   247017  35403   198744  107562    175602    130649    175895];
names3 = {'',    '',    '',   'ARRDC3','TP63','HOMER1', '',     'VCL',  '',     'CXCL12', 'CCDC85B','CYP2E1','PLEKHF2'};
allids{3} = ids3;
allnames{3} = names3;



%---------------------------------------------------
% Plot the series
%---------------------------------------------------


f=fopen('biologistplots.tex','w');

%-------------------------------------------
% write a LaTeX file that will display the plots
%-------------------------------------------
  fprintf(f,'\\documentclass{article}\n')
  fprintf(f,'\\usepackage{graphicx}\n')
  fprintf(f,'\\addtolength{\\textwidth}{5cm}\n');
  fprintf(f,'\\addtolength{\\hoffset}{-2.0cm}\n');
  fprintf(f,'\\begin{document}\n');



timeshift=300;

for indtemp=1:3,

ids=allids{indtemp};
names=allnames{indtemp};
maxk=length(ids);
k=1;while k<=maxk,
  k
  l = find(allresults_ensemblids == ids(k));
  model=allresults_jointmodels{l};
  if ~isempty(model),
    % The next three lines are needed to fix long variable names
    % which Octave does not save correctly.
    model.disimdelaytransformationsettings=model.disimdelaytransformationsetting;
    model.disimvariancetransformationsettings=model.disimvariancetransformationsett;
    model.disimdecaytransformationsettings=model.disimdecaytransformationsetting;

    ensemblid=allresults_ensemblids(l);
    loglik=allresults_loglikelihoods(l,3);
    forlatex=1;
    plottitle=makeplottitle(model,loglik,ensemblid,names{k},forlatex);
  
    predicttimes=timeshift + 1280*(([0:256]'/256).^2);

if 1,
    plotpredictions(model,predicttimes,2,1,0,1,1,0,[],min(model.t{1}));
% (gpsim3model,predicttimes,timescaletype,stdevmultiplier,plotrbf,plotpol,plotrna,plotpolrna,plottitle,timeshift);  

    if 0,
    h=axes;
    ylim=get(h,'ylim');
    xlim=get(h,'xlim');
    xloc=xlim(1);
    %yloc=ylim(1)+(ylim(2)-ylim(1))*2.6238;
    %yloc=ylim(1)+(ylim(2)-ylim(1))*2.4238;
    yloc=ylim(1)+(ylim(2)-ylim(1))*1.4238;
    h=text(xloc,yloc,plottitle);
    set(h,'fontname','Helvetica');
    set(h,'fontsize',20);
    axis off;
    drawnow;
    end;  

    % Octave version
    % print(sprintf('ENSG%011d_GP_jointfit.eps'),'-depsc','-S700,600');

    % Matlab version
    h=gcf;
    set(h,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
    print(sprintf('ENSG%011d_GP_jointfit.eps',ensemblid),'-depsc','-r400');
end;

    if ((k>1) || (indtemp>1)), fprintf(f,'\\newpage\n'); end;
    %fprintf(f,'$\\phantom{a}$\n');
    %fprintf(f,'\\hspace{-2cm}\n');
    fprintf(f,'\\noindent %s',plottitle);
    fprintf(f,'\\vspace{1cm}\n');
    fprintf(f,'\n\n');
    fprintf(f,'\\includegraphics[width=0.9\\textwidth]{%s}\n', sprintf('ENSG%011d_GP_jointfit.eps',ensemblid)  );    
    
  end;

  k=k+1; 
end;

end;


  fprintf(f,'\\end{document}');
  fclose(f);
