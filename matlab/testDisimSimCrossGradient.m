function [diffgradient,kernelgradient]=testDisimSimCrossGradient(disimkernel,simkernel,t,covGrad);

if isempty(covGrad),
  covGrad=randn(2*kernel.numBlocks);
  covGrad=covGrad+covGrad';
end;

[disimpars,nams]=kernExtractParam(disimkernel);
disimkernel=kernExpandParam(disimkernel,disimpars);

[simpars,nams]=kernExtractParam(simkernel);
simkernel=kernExpandParam(simkernel,simpars);

pars=[disimpars(1) disimpars(2) disimpars(3) disimpars(4) disimpars(5) disimpars(6) disimpars(7)];

kernelgradient=disimXsimKernGradient(disimkernel,simkernel,t,covGrad);

diffgradient=zeros(1,length(pars));
for k=1:length(pars),
  parstemp=pars;
  
  if (abs(parstemp(k))>1e-6)
    myeps=1e-6;
  else
    myeps=abs(parstemp(k))*0.5;
  end;

  parstemp(k)=pars(k)-myeps;
  disimkerntemp=kernExpandParam(disimkernel,parstemp);
  simkerntemp=kernExpandParam(simkernel,parstemp(1:3));  
  Ktemp1=disimXsimKernCompute(disimkerntemp,simkerntemp,t);
  
  parstemp(k)=pars(k)+myeps;
  disimkerntemp=kernExpandParam(disimkernel,parstemp);
  simkerntemp=kernExpandParam(simkernel,parstemp(1:3));  
  Ktemp2=disimXsimKernCompute(disimkerntemp,simkerntemp,t);
  
  diffgradient(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));
end;

[kernelgradient' diffgradient']
