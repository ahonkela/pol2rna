function plot_model(model, ensemblid, loglik, symbol),

% symbols = r.symbols(A);
% plot_model(results_jointmodels{k}, results_ensemblids(k), results_loglikelihoods(k), symbols{k})

timeshift = 300;

    
%k=2;
%model = results_jointmodels{k};
%ensemblid=results_ensemblids(k);
%loglik=results_loglikelihoods(k);
plottitle=makeplottitle(model,loglik,ensemblid,symbol,0);
predicttimes=timeshift + 1280*(([0:256]'/256).^2);
plotpredictions(model,predicttimes,2,1,0,1,1,0,[],min(model.t{1}));
h=gca;
ylim=get(h,'ylim');
xlim=get(h,'xlim');
%xloc=xlim(1);
%yloc=ylim(1)+(ylim(2)-ylim(1))*2.6238;
%h=text(xloc,yloc,plottitle);
subplot(2, 1, 1);
title(plottitle)
set(h,'fontname','Helvetica');
%drawnow;
%print(sprintf('ENSG%011d_joint_GP_2012-05-07.eps',ensemblid),'-depsc','-S1024,800');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [1024/72, 786/72]);
print(sprintf('ENSG%011d_joint_GP_2012-05-07.png',ensemblid),'-dpng','-r72');
