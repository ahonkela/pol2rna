function [difftotalgradient,modeltotalgradient]=testModelGradientBasic(model,updatekernel,updatemean);


[pars,nams]=gpdisimExtractParam(model);
model=gpdisimExpandParam(model,pars);

modeltotalgradient=gpdisimLogLikeGradients(model,updatekernel,updatemean);

myepses=[1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];
myeps=1e-2;
for k=1:length(pars),
  k
  bestgraddiff=inf;
  
  for l=1:length(myepses),
    parstemp=pars;
    
    parstemp(k)=pars(k)-myeps;
    modeltemp=gpdisimExpandParam(model,parstemp,updatekernel,updatemean);
    valtemp1=gpdisimLogLikelihood(modeltemp);
    
    parstemp(k)=pars(k)+myeps;
    modeltemp=gpdisimExpandParam(model,parstemp,updatekernel,updatemean);
    valtemp2=gpdisimLogLikelihood(modeltemp);

    tempdiffgradient=(valtemp2-valtemp1)/(2*myeps);
    
    tempgraddiff=abs(modeltotalgradient(k)-tempdiffgradient);
    if (tempgraddiff<bestgraddiff),
      bestgraddiff=tempgraddiff;
      difftotalgradient(k)=tempdiffgradient;
    end;
  end;
end;

[modeltotalgradient' difftotalgradient']


%dim=length(model.t);
%K=model.K;
%varpar=model.kern.comp{1}.comp{1}.variance;
%(-0.5*dim/varpar + 0.5/varpar*model.m'*inv(model.K)*model.m)*sigmoidabTransform(varpar,'gradfact',[0 10])