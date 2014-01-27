function [pol2predictions,pol2vars,rnapredictions,rnavars]=plotpredtemp(gpsim3model,predicttimes,timescaletype,stdevmultiplier,plotdependent,plottitle);

fontsize=8;
mylinewidth=2;

if isempty(predicttimes),
  predicttimes=[0:1:1280]';
end;

numGenes=gpsim3model.numGenes;
%if plotrna==0,
%  numGenes=0;
%end;


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
%  kwidth=length(predicttimes);
  kwidth=50;
  temptimes=predicttimes(k:min([length(predicttimes) k+kwidth-1]));
  [temppriormeans,tempmeans,covmatrix,rbfmeans,rbfcovmatrix]=gpasimTemp4Predict(gpsim3model,temptimes,plotdependent);
  tempdiag=diag(covmatrix);

  plength=length(temptimes);
  
  pol2predictions(k:k+plength-1)=tempmeans(1:plength);
  pol2vars(k:k+plength-1)=tempdiag(1:plength);

  rbfdiag=diag(rbfcovmatrix);
  rbfpredictions(k:k+plength-1)=rbfmeans(1:plength);
  rbfvars(k:k+plength-1)=rbfdiag(1:plength);

  if plotdependent==1,
    rnapredictions(k:k+plength-1)=tempmeans(plength+1:2*plength);
    rnavars(k:k+plength-1)=tempdiag(plength+1:2*plength);
  else
    [temppriormeans,tempmeans,covmatrix] = gpasimPredictRNAOnly(gpsim3model,temptimes);
    rnapredictions(k:k+plength-1)=tempmeans(1:plength);
    tempdiag=diag(covmatrix);
    rnavars(k:k+plength-1)=tempdiag(1:plength);
  end;

  if (plotdependent==1) | (plotdependent==2),
    [temppriormeans,tempmeans,covmatrix]=gpasimPredictRNAConditional(gpsim3model,temptimes);
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




%pol2predictions=posteriormeans(1:length(predicttimes));
%pol2vars=covdiag(1:length(predicttimes));

%rnapredictions=posteriormeans(length(predicttimes)+1:2*length(predicttimes));
%rnavars=covdiag(length(predicttimes)+1:2*length(predicttimes));

%for k=1:length(predicttimes),
%  k
%  [predmeans,predcov]=gpasimTemp3Predict(gpsim3model,predicttimes(k),predicttimes(k));
%  pol2predictions(k)=predmeans(1);
%  rnapredictions(k)=predmeans(2);
%  pol2vars(k)=predcov(1,1);
%  rnavars(k)=predcov(2,2);
%end;



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
end;
  
% I=find(rbfvars<0);rbfvars(I)=nan;
% I=find(pol2vars<0);pol2vars(I)=nan;
% I=find(rnavars<0);rnavars(I)=nan;


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
  ticklabels{k}=sprintf('%d',ticktimes_notransformation(k));
end;



clf;


%--------------------------------------
% Plot RNA predicted series
%--------------------------------------
if plotdependent==1,
  subplot(4,1,4);
end;
  
if ((plotdependent==1) | (plotdependent==2)),
  
  rnavals=gpsim3model.y(length(pol2times)+1:length(pol2times)+length(rnatimes));
  
  % hold on; plot(predicttimes,rnapredictions+stdevmultiplier*(rnavars.^0.5),'c-');
  % hold on; plot(predicttimes,rnapredictions-stdevmultiplier*(rnavars.^0.5),'m-');
  
  hold on;
  tempc=[0.5 1.0 0.9];  
  tempy=rnapredictionsconditional;
  tempd=stdevmultiplier*(rnavarsconditional.^0.5);

  % tempy
  % tempd
  % pause
  
  % pause
  tempd=real(tempd);
  % pause
  
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
  
end;




%figure;
%subplot(2,1,1);
%plot(cumsum(rbfpredictions),'k-');
%subplot(2,1,2);
%plot(pol2predictions,'b-');
