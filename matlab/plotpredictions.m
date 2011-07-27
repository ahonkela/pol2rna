function [pol2predictions,pol2vars,rnapredictions,rnavars]=plotpredictions(gpsim3model,predicttimes,timescaletype,stdevmultiplier,plotdependent,plottitle);

fontsize=9;
mylinewidth=2;

if isempty(predicttimes),
  predicttimes=[0:1:1280]';
end;

numGenes=gpsim3model.numGenes;


% predict only a few points at a time to avoid numerical problems
pol2predictions=nan*predicttimes;
pol2vars=nan*predicttimes;
rnapredictions=nan*predicttimes;
rnavars=nan*predicttimes;
rnapredictionsconditional=nan*predicttimes;
rnavarsconditional=nan*predicttimes;
rbfpredictions=nan*predicttimes;
rbfvars=nan*predicttimes;
covdiag=nan*[predicttimes;predicttimes];
k=1;
while k<length(predicttimes),
  k
  kwidth=length(predicttimes);
  %kwidth=50;
  temptimes=predicttimes(k:min([length(predicttimes) k+kwidth-1]));
  [temppriormeans,tempmeans,covmatrix,rbfmeans,rbfcovmatrix]=gpnddisimPredict(gpsim3model,temptimes,plotdependent);
  tempdiag=diag(covmatrix);

  plength=length(temptimes);
  
  pol2predictions(k:k+plength-1)=tempmeans(1:plength);
  pol2vars(k:k+plength-1)=tempdiag(1:plength);

  if ~isempty(rbfmeans),
    rbfdiag=diag(rbfcovmatrix);
    rbfpredictions(k:k+plength-1)=rbfmeans(1:plength);
    rbfvars(k:k+plength-1)=rbfdiag(1:plength);
  end;
    
  if plotdependent==1,
    rnapredictions(k:k+plength-1)=tempmeans(plength+1:2*plength);
    rnavars(k:k+plength-1)=tempdiag(plength+1:2*plength);
  else
    [temppriormeans,tempmeans,covmatrix] = gpnddisimPredictRNAOnly(gpsim3model,temptimes);
    rnapredictions(k:k+plength-1)=tempmeans(1:plength);
    tempdiag=diag(covmatrix);
    rnavars(k:k+plength-1)=tempdiag(1:plength);
  end;

  if ((plotdependent==1) | (plotdependent==2) | (plotdependent==3)),
    [temppriormeans,tempmeans,covmatrix]=gpnddisimPredictRNAConditional(gpsim3model,temptimes);
    tempdiag=diag(covmatrix);
    
    rnapredictionsconditional(k:k+plength-1)=tempmeans(1:plength);
    rnavarsconditional(k:k+plength-1)=tempdiag(1:plength);
  end;
  
  k=k+kwidth;
end;

rbfvars=real(rbfvars);
I=find(rbfvars<0);
rbfvars(I)=0;

pol2vars=real(pol2vars);
I=find(pol2vars<0);
pol2vars(I)=0;

rnavars=real(rnavars);
I=find(rnavars<0);
rnavars(I)=0;

rnavarsconditional=real(rnavarsconditional);
I=find(rnavarsconditional<0);
rnavarsconditional(I)=0;


if 0,
  rbfpredictions
  pause
  pol2predictions
  pause
  rnapredictions
  pause
  rbfvars
  pause
  pol2vars
  pause
  rnavars
  pause
  rnavarsconditional
  pause

  I=find(rbfvars<0);rbfvars(I)=nan;
  I=find(pol2vars<0);pol2vars(I)=nan;
  I=find(rnavars<0);rnavars(I)=nan;
end;
  


ticktimes_notransformation=gpsim3model.t;
if timescaletype==0,
  temptimes=predicttimes;
  pol2times=gpsim3model.t;
  rnatimes=gpsim3model.t;
  ticktimes=ticktimes_notransformation;
elseif timescaletype==1,
  temptimes=log(predicttimes+1);
  pol2times=log(gpsim3model.t+1);
  rnatimes=log(gpsim3model.t+1);
  ticktimes=log(ticktimes_notransformation+1);
elseif timescaletype==2,
  temptimes=sqrt(predicttimes);
  pol2times=sqrt(gpsim3model.t);
  rnatimes=sqrt(gpsim3model.t);
  ticktimes=sqrt(ticktimes_notransformation);
end;
ticklabels={};
for k=1:length(ticktimes), 
  ticklabels{k}=sprintf('%.2f',ticktimes_notransformation(k));
end;



clf;

%--------------------------------------
% Plot Driving RBF
%--------------------------------------
drawme=0;
showtitle=0;
if plotdependent==1,
  subplot(4,1,1);
  drawme=1;
  showtitle=1;
elseif plotdependent==0
  subplot(3,1,1);
  drawme=1;
  showtitle=1;
end;
drawme=0;

if drawme==1,
  
  hold on;
  tempc=[0.8 0.8 0.8];  
  tempy=rbfpredictions;
  tempd=stdevmultiplier*(real(rbfvars.^0.5));
  % tempy
  % tempd
  
  try
    h=fill([temptimes; temptimes(end:-1:1)], [tempy+tempd; tempy(end:-1:1)-tempd(end:-1:1)],tempc);
    set(h,'EdgeColor',tempc);
    set(h,'FaceAlpha',0.5);
    set(h,'EdgeAlpha',0.5);
  catch
    hold on; h=plot(temptimes,tempy+tempd,'k--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
    hold on; h=plot(temptimes,tempy-tempd,'k--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
  end;
  
  
  hold on; h=plot(temptimes,rbfpredictions,'k-');set(h,'LineWidth',mylinewidth);
  
  h=gca;
  set(h,'xtick',ticktimes);
  set(h,'xticklabel',ticklabels);
  set(h,'fontsize',fontsize);
  
  ylabel('driving RBF');
  
  if showtitle==1,
    if ~isempty(plottitle),
      title(plottitle);
    end;  
  end;
  
end;


%--------------------------------------
% Plot Input series
%--------------------------------------
drawme=0;
showtitle=1;
if plotdependent==1,
  subplot(3,1,1);
  drawme=1;
elseif plotdependent==0,
  subplot(3,1,2);
  drawme=1;
elseif plotdependent==3,
  subplot(3,1,1);
  drawme=1;
  showtitle=1;
end;

if drawme==1,
  
  pol2vals=gpsim3model.y(1:length(pol2times));
  
  hold on;
  tempc=[1.0 0.7 0.7];  
  tempy=pol2predictions;
  tempd=stdevmultiplier*(pol2vars.^0.5);
  
  try
    h=fill([temptimes; temptimes(end:-1:1)], [tempy+tempd; tempy(end:-1:1)-tempd(end:-1:1)],tempc);
    set(h,'EdgeColor',tempc);
    set(h,'FaceAlpha',0.5);
    set(h,'EdgeAlpha',0.5);    
  catch  
    hold on; h=plot(temptimes,tempy+tempd,'r--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
    hold on; h=plot(temptimes,tempy-tempd,'r--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
  end;
  
  hold on; h=plot(temptimes,pol2predictions,'r-');set(h,'LineWidth',mylinewidth);
  
  hold on; h=plot(pol2times,pol2vals,'k:o');
    
  h=gca;
  set(h,'xtick',ticktimes);
  set(h,'xticklabel',ticklabels);
  set(h,'fontsize',fontsize);
  
  ylabel('input-series');

  if showtitle==1,
    if ~isempty(plottitle),
      title(plottitle);
    end;  
  end;
  
end;



%--------------------------------------
% Plot RNA series
%--------------------------------------
drawme=0;
showtitle=0;
if plotdependent==1,
  subplot(3,1,2);
  drawme=1;
elseif plotdependent==0,
  subplot(3,1,3);
  drawme=1;
elseif plotdependent==3,
  subplot(3,1,2);
  drawme=1;
end;

if drawme==1,
  
  rnavals=gpsim3model.y(length(pol2times)+1:length(pol2times)+length(rnatimes));
  
  % hold on; plot(predicttimes,rnapredictions+stdevmultiplier*(rnavars.^0.5),'c-');
  % hold on; plot(predicttimes,rnapredictions-stdevmultiplier*(rnavars.^0.5),'m-');
  
  hold on;
  tempc=[0.7 1.0 0.7];  
  tempy=rnapredictions;
  tempd=stdevmultiplier*(rnavars.^0.5);
  
  try
    h=fill([temptimes; temptimes(end:-1:1)], [tempy+tempd; tempy(end:-1:1)-tempd(end:-1:1)],tempc);
    set(h,'EdgeColor',tempc);
    set(h,'FaceAlpha',0.5);
    set(h,'EdgeAlpha',0.5);  
  catch
    hold on; h=plot(temptimes,tempy+tempd,'g--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
    hold on; h=plot(temptimes,tempy-tempd,'g--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
  end;
  
  
  % rnapredictions
  % pause
  
  hold on; h=plot(temptimes,rnapredictions,'g-');set(h,'LineWidth',mylinewidth);
  
  hold on; plot(rnatimes,rnavals,'k:o');
  
  %axis([min(predicttimes) max(predicttimes) min(rnavals)-sqrt(var(rnavals)) max(rnavals)+sqrt(var(rnavals))]);
  
  h=gca;
  set(h,'xtick',ticktimes);
  set(h,'xticklabel',ticklabels);
  set(h,'fontsize',fontsize);
  
  ylabel('RNA-series');  
  
  if showtitle==1,
    if ~isempty(plottitle),
      title(plottitle);
    end;  
  end;
  
end;


%--------------------------------------
% Plot RNA predicted series
%--------------------------------------
drawme=0;
showtitle=0;
if plotdependent==1,
  subplot(3,1,3);
  drawme=1;
elseif plotdependent==2,
  drawme=1;
  showtitle=1;
elseif plotdependent==3,
  subplot(3,1,3);
  drawme=1;
end;
  
if (drawme==1),
  
  rnavals=gpsim3model.y(length(pol2times)+1:length(pol2times)+length(rnatimes));
  
  % hold on; plot(predicttimes,rnapredictions+stdevmultiplier*(rnavars.^0.5),'c-');
  % hold on; plot(predicttimes,rnapredictions-stdevmultiplier*(rnavars.^0.5),'m-');
  
  hold on;
  tempc=[0.5 1.0 0.9];  
  tempy=rnapredictionsconditional;
  tempd=stdevmultiplier*(rnavarsconditional.^0.5);
  
  
  % pause
  tempd=real(tempd)
  %pause
  
  try
    h=fill([temptimes; temptimes(end:-1:1)], [tempy+tempd; tempy(end:-1:1)-tempd(end:-1:1)],tempc);
    set(h,'EdgeColor',tempc);
    set(h,'FaceAlpha',0.5);
    set(h,'EdgeAlpha',0.5);  
  catch
    hold on; h=plot(temptimes,tempy+tempd,'c--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
    hold on; h=plot(temptimes,tempy-tempd,'c--');set(h,'Color',tempc);set(h,'LineWidth',mylinewidth);
  end;
  
  
  hold on; h=plot(temptimes,rnapredictionsconditional,'c-');set(h,'LineWidth',mylinewidth);
  
  hold on; plot(rnatimes,rnavals,'k:o');
  
  %axis([min(predicttimes) max(predicttimes) min(rnavals)-sqrt(var(rnavals)) max(rnavals)+sqrt(var(rnavals))]);
  
  h=gca;
  set(h,'xtick',ticktimes);
  set(h,'xticklabel',ticklabels);
  set(h,'fontsize',fontsize);
  
  ylabel('RNA predicted from input');  
  
  if showtitle==1,
    if ~isempty(plottitle),
      title(plottitle);
    end;  
  end;
  
end;




%figure;
%subplot(2,1,1);
%plot(cumsum(rbfpredictions),'k-');
%subplot(2,1,2);
%plot(pol2predictions,'b-');
