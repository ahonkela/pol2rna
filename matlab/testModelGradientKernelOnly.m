function [difftotalgradient,modeltotalgradient]=testModelGradientKernelOnly(model,updatekernel,updatemean);


[pars,nams]=gpdisimExtractParam(model);
model=gpdisimExpandParam(model,pars);
dim = size(model.y, 1);
m=model.m;
K=model.K;
invK=model.invK;
logDetK=model.logDetK;

updatekernel=1;
updatemean=0;
modeltotalgradient=gpdisimLogLikeGradients(model,updatekernel,updatemean);

myeps=1e-5;
for k=1:length(pars),
  parstemp=pars;
  
  parstemp(k)=pars(k)-myeps;
  modeltemp=gpdisimExpandParam(model,parstemp,updatekernel,updatemean);
  valtemp1=gpdisimLogLikelihood(modeltemp);
  
  parstemp(k)=pars(k)+myeps;
  modeltemp=gpdisimExpandParam(model,parstemp,updatekernel,updatemean);
  valtemp2=gpdisimLogLikelihood(modeltemp);
  
  difftotalgradient(k)=(valtemp2-valtemp1)/(2*myeps);


  parstemp(k)=pars(k)-myeps;
  modeltemp=gpdisimExpandParam(model,parstemp,1,1);
  valtemp1=0.5*(-dim*log(2*pi) - modeltemp.logDetK - m'*modeltemp.invK*m);
  
  parstemp(k)=pars(k)+myeps;
  modeltemp=gpdisimExpandParam(model,parstemp,1,1);
  valtemp2=0.5*(-dim*log(2*pi) - modeltemp.logDetK - m'*modeltemp.invK*m);

  difftotalgradient2(k)=(valtemp2-valtemp1)/(2*myeps);
end;

modeltotalgradient
difftotalgradient
difftotalgradient2


%dim=length(model.t);
%K=model.K;
%varpar=model.kern.comp{1}.comp{1}.variance;
%(-0.5*dim/varpar + 0.5/varpar*model.m'*inv(model.K)*model.m)*sigmoidabTransform(varpar,'gradfact',[0 10])