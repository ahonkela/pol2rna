function plotprofiles(ensemblid,tempprofile,firstk,lastk,centerofmass1);
clf;
timepoints=[0 5 10 20 40 80 160 320 640 1280];
ymax=max(max(tempprofile));
for k=firstk:lastk,
  subplot(lastk-firstk+1,1,k-firstk+1);
  h=bar([1:size(tempprofile,2)]*200,tempprofile(k,:));
  set(h(1),'barwidth',1.0)
  h=gca;
  set(h,'ylim',[0 ymax]);
  if k<lastk, 
    set(h,'xtick',[]);
  end;
  if k==1,
    title(sprintf('ENSEMBL id %d',ensemblid));
  end;

  if ~isempty(centerofmass1),
    hold on;
    h=plot([centerofmass1(k) centerofmass1(k)],[0 ymax]);
    set(h(1),'linewidth',3.0);
    set(h(1),'color',[1 0 0]);
  end;
ylabel(sprintf('%dmin',timepoints(k)));
end;

