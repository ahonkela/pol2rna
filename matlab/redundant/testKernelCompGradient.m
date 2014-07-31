function [diffgradient,kernelgradient]=testKernelCompGradient(kernel,t,covGrad,compindex);

if isempty(covGrad),
  covGrad=randn(length(t)*kernel.numBlocks);
  covGrad=covGrad+covGrad';
end;

[pars,nams]=kernExtractParam(kernel);
kernel=kernExpandParam(kernel,pars);
[pars,nams]=kernExtractParam(kernel)
pause


K=kernCompute(kernel,t);
kernelgradient=kernGradient(kernel,t,covGrad);

diffgradient=zeros(1,length(pars));
for k=1:length(pars),
  parstemp=pars;
  
  if (abs(parstemp(k))>1e-6)
    myeps=1e-6;
  else
    myeps=abs(parstemp(k))*0.5;
  end;

  parstemp(k)=pars(k)-myeps;
  kerntemp=kernExpandParam(kernel,parstemp);
  %kerntemp
  %pause  
  Ktemp1=kernCompute(kerntemp.comp{compindex},t);
  
  parstemp(k)=pars(k)+myeps;
  kerntemp=kernExpandParam(kernel,parstemp);
  Ktemp2=kernCompute(kerntemp.comp{compindex},t);
  
  diffgradient(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));
end;

kernelgradient 
diffgradient
