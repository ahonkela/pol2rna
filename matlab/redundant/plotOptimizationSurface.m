function lls=plotOptimizationSurface(model,transforminfo,paramind1,paramind2);


param_x1=transforminfo.settings{paramind1}(1);
param_x2=transforminfo.settings{paramind1}(2);
param_y1=transforminfo.settings{paramind2}(1);
param_y2=transforminfo.settings{paramind2}(2);

[pars,nams]=gpdisimExtractParam(model);
parst=transformParametersWithSettings(pars,transforminfo,'atox');

npointsx=10;
npointsy=10;

lls=zeros(npointsy,npointsx);
for i=1:npointsx,
  i
  for j=1:npointsy,
    parsttemp=parst;
    parsttemp(paramind1)=param_x1+(i-1)*(param_x2-param_x1)/(npointsx-1);
    parsttemp(paramind2)=param_y1+(j-1)*(param_y2-param_y1)/(npointsy-1);
    
    parstemp=transformParametersWithSettings(parsttemp,transforminfo,'xtoa');
    modeltemp=gpdisimExpandParam(model,parstemp);
    lltemp=gpdisimLogLikelihood(modeltemp);
    lls(j,i)=lltemp;
  end;
end;


llpower=0.5;
I=find(lls(:)>0);
lls(I)=lls(I).^llpower;
I=find(lls(:)<0);
lls(I)=-(-lls(I)).^llpower;


imagesc(lls);
xlabel(nams{paramind1});
ylabel(nams{paramind2});

xticknames=cell(npointsx,1);
for i=1:npointsx,
  xval=param_x1+(i-1)*(param_x2-param_x1)/(npointsx-1);
  xticknames{i}=sprintf('%e',xval);
end;
yticknames=cell(npointsy,1);
for i=1:npointsy,
  yval=param_y1+(i-1)*(param_y2-param_y1)/(npointsy-1);
  yticknames{i}=sprintf('%e',yval);
end;

h=gca;
set(h,'xticklabel',xticknames);
set(h,'yticklabel',yticknames);

