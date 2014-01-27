function [difftotalgradient,modeltotalgradient]=testModelGradient(model,updatekernel,updatemean);

[pars,nams]=gpdisimExtractParam(model);
model=gpdisimExpandParam(model,pars);



if 0,

disimkern=model.kern.comp{1}.comp{2}; 
covGrad=ones(length(model.t));
modelkerngradient1=disimKernGradient(disimkern,model.t,covGrad);

%[pars,nams]=kernExtractParam(disimkern);
pars=[disimkern.di_decay disimkern.inverseWidth disimkern.di_variance disimkern.decay disimkern.variance disimkern.rbf_variance];

for k=1:length(pars),
  parstemp=pars;
  myeps=1e-5;
  
  parstemp(k)=pars(k)-myeps;
  kerntemp=disimKernExpandParam(disimkern,parstemp);
  Ktemp1=kernCompute(kerntemp,model.t);
    
  parstemp(k)=pars(k)+myeps;
  kerntemp=disimKernExpandParam(disimkern,parstemp);
  Ktemp2=kernCompute(kerntemp,model.t);

  diffkerngradient1(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));
end;
modelkerngradient1
diffkerngradient1

pause

end;


if 0,

%covGrad=ones(length(model.t));
%covGrad = -model.invK(10:18,10:18) + model.invK(10:18,10:18)*model.m(10:18)*model.m(10:18)'*model.invK(10:18,10:18);
covGrad=randn(9);
disimkern=model.kern.comp{1}.comp{2}; 
modelkerngradient2=kernGradient(disimkern,model.t,covGrad);
[pars,nams]=kernExtractParam(disimkern);

for k=1:length(pars),
  parstemp=pars;
  myeps=1e-5;
  
  parstemp(k)=pars(k)-myeps;
  kerntemp=kernExpandParam(disimkern,parstemp);
  Ktemp1=kernCompute(kerntemp,model.t);
    
  parstemp(k)=pars(k)+myeps;
  kerntemp=kernExpandParam(disimkern,parstemp);
  Ktemp2=kernCompute(kerntemp,model.t);

  diffkerngradient2(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));
end;
modelkerngradient2
diffkerngradient2

pause

end;



if 0,
  
covGrad = randn(length(model.t));
simkern=model.kern.comp{1}.comp{1}; 
disimkern=model.kern.comp{1}.comp{2}; 
modelkerngradient2a=disimXsimKernGradient(disimkern,simkern,model.t,covGrad);
pars=[disimkern.di_decay disimkern.inverseWidth disimkern.di_variance disimkern.decay disimkern.variance disimkern.rbf_variance];

for k=1:length(pars),
  parstemp=pars;
  myeps=1e-5;

  parstemp(k)=pars(k)-myeps;
  disimkerntemp=disimKernExpandParam(disimkern,[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  simkerntemp=simKernExpandParam(simkern,[parstemp(1) parstemp(2) parstemp(3)]);  
  Ktemp1=disimXsimKernCompute(disimkerntemp,simkerntemp,model.t);
    
  parstemp(k)=pars(k)+myeps;
  disimkerntemp=disimKernExpandParam(disimkern,[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  simkerntemp=simKernExpandParam(simkern,[parstemp(1) parstemp(2) parstemp(3)]);
  Ktemp2=disimXsimKernCompute(disimkerntemp,simkerntemp,model.t);

  diffkerngradient2a(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));
end;
modelkerngradient2a
diffkerngradient2a

pause

end;


effectkern=model.kern.comp{1};
simkern=model.kern.comp{1}.comp{1}; 
disimkern=model.kern.comp{1}.comp{2};
noisekern=model.kern.comp{2};
whitekern1=model.kern.comp{2}.comp{1};
whitekern2=model.kern.comp{2}.comp{2};

randn('seed',1252823532);
covGrad = randn(length(model.t)*2);
covGrad = covGrad+covGrad';
%modelkerngradient2a2=kernGradient(model.kern,model.t,covGrad);

%trace(covGrad(1:9,1:9))
%trace(covGrad(10:18,10:18))

%covGrad = randn(length(model.t)*1);
%modelkerngradient2a2=kernGradient(model.kern.comp{2}.comp{2},model.t,covGrad(10:18,10:18));
%modelkerngradient2a2=disimXsimKernGradient(model.kern.comp{1}.comp{2},model.kern.comp{1}.comp{1},model.t,covGrad);

pars=[disimkern.di_decay disimkern.inverseWidth disimkern.di_variance disimkern.decay disimkern.variance disimkern.rbf_variance whitekern1.variance whitekern2.variance]

effectK=kernCompute(effectkern,model.t);
sum(sum(effectK(1:9,1:9).*covGrad(1:9,1:9)))

modelkerngradient2a2=kernGradient(model.kern.comp{1},model.t,covGrad)

%whiteKernGradient(model.kern.comp{2}.comp{1},model.t,covGrad(1:9,1:9))
%whiteKernGradient(model.kern.comp{2}.comp{2},model.t,covGrad(10:18,10:18))

diffkerngradient2a2=zeros(1,length(pars));
for k=3:3, %1:length(pars),
  parstemp=pars;
  myeps=1e-6;

  parstemp(k)=pars(k);
  disimkerntemp=disimKernExpandParam(disimkern,[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  simkerntemp=simKernExpandParam(simkern,[parstemp(1) parstemp(2) parstemp(3)]);  
  whitekern1temp=whiteKernExpandParam(whitekern1,[parstemp(7)]);
  whitekern2temp=whiteKernExpandParam(whitekern2,[parstemp(8)]);
  Ktemp1a=simKernCompute(simkerntemp,model.t);
  Ktemp1a2=whiteKernCompute(whitekern1temp,model.t);
  Ktemp1b=disimXsimKernCompute(disimkerntemp,simkerntemp,model.t);
  Ktemp1c=disimKernCompute(disimkerntemp,model.t);
  Ktemp1c2=whiteKernCompute(whitekern2temp,model.t);
  Ktemp1=[Ktemp1a Ktemp1b';Ktemp1b Ktemp1c];
  fprintf(1,'Simple computation of di_variance gradient\n')
%  sum(sum((Ktemp1/disimkern.di_variance).*covGrad))
  sum(sum((Ktemp1/disimkern.di_variance).*covGrad))*sigmoidabTransform(pars(3),'gradfact',disimkern.transforms(3).transformsettings)
  
  
  parstemp(k)=pars(k)-myeps;
  disimkerntemp=disimKernExpandParam(disimkern,[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  simkerntemp=simKernExpandParam(simkern,[parstemp(1) parstemp(2) parstemp(3)]);  
  whitekern1temp=whiteKernExpandParam(whitekern1,[parstemp(7)]);
  whitekern2temp=whiteKernExpandParam(whitekern2,[parstemp(8)]);
  Ktemp1a=simKernCompute(simkerntemp,model.t);
  Ktemp1a2=whiteKernCompute(whitekern1temp,model.t);
  Ktemp1b=disimXsimKernCompute(disimkerntemp,simkerntemp,model.t);
  Ktemp1c=disimKernCompute(disimkerntemp,model.t);
  Ktemp1c2=whiteKernCompute(whitekern2temp,model.t);
  Ktemp1=[Ktemp1a Ktemp1b';Ktemp1b Ktemp1c];
  %Ktemp1=[Ktemp1a2 zeros(length(model.t));zeros(length(model.t)) Ktemp1c2];
  %Ktemp1=full(Ktemp1c2);

%  tempkern=kernCompute(model.kern.comp{1},model.t);
%  tempkern-Ktemp1
  
  parstemp(k)=pars(k)+myeps;
  disimkerntemp=disimKernExpandParam(disimkern,[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  simkerntemp=simKernExpandParam(simkern,[parstemp(1) parstemp(2) parstemp(3)]);  
  whitekern1temp=whiteKernExpandParam(whitekern1,[parstemp(7)]);
  whitekern2temp=whiteKernExpandParam(whitekern2,[parstemp(8)]);
  Ktemp2a=simKernCompute(simkerntemp,model.t);
  Ktemp2a2=whiteKernCompute(whitekern1temp,model.t);
  Ktemp2b=disimXsimKernCompute(disimkerntemp,simkerntemp,model.t);
  Ktemp2c=disimKernCompute(disimkerntemp,model.t);
  Ktemp2c2=whiteKernCompute(whitekern2temp,model.t);
  Ktemp2=[Ktemp2a Ktemp2b';Ktemp2b Ktemp2c];
  %Ktemp2=[Ktemp2a2 zeros(length(model.t));zeros(length(model.t)) Ktemp2c2];
  %Ktemp2=full(Ktemp2c2);  
  
  diffkerngradient2a2(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));
end;

fprintf('diffkerngradient2a2 before factors:\n')
diffkerngradient2a2

factors(1)=sigmoidabTransform(pars(1),'gradfact',disimkern.transforms(1).transformsettings);
factors(2)=sigmoidabTransform(pars(2),'gradfact',disimkern.transforms(2).transformsettings);
factors(3)=sigmoidabTransform(pars(3),'gradfact',disimkern.transforms(3).transformsettings);
factors(4)=sigmoidabTransform(pars(4),'gradfact',disimkern.transforms(4).transformsettings);
factors(5)=sigmoidabTransform(pars(5),'gradfact',disimkern.transforms(5).transformsettings);
factors(6)=sigmoidabTransform(pars(6),'gradfact',disimkern.transforms(6).transformsettings);
factors(7)=sigmoidabTransform(pars(7),'gradfact',whitekern1.transforms(1).transformsettings);
factors(8)=sigmoidabTransform(pars(8),'gradfact',whitekern2.transforms(1).transformsettings);
diffkerngradient2a2=diffkerngradient2a2.*factors;

modelkerngradient2a2
diffkerngradient2a2

pause





if 0,

covGrad = randn(length(model.t));
simkern=model.kern.comp{1}.comp{1}; 
disimkern=model.kern.comp{1}.comp{2}; 
modelkerngradient2b=disimXsimKernGradient(disimkern,simkern,model.t,covGrad);
[pars,nams]=kernExtractParam(model.kern);

for k=1:length(pars),
  parstemp=pars;
  myeps=1e-5;

  tempkern=model.kern;
  parstemp(k)=pars(k)-myeps;
  disimkerntemp=kernExpandParam(disimkern,[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  simkerntemp=kernExpandParam(simkern,[parstemp(1) parstemp(2) parstemp(3)]);  
  Ktemp1=disimXsimKernCompute(disimkerntemp,simkerntemp,model.t);
    
  parstemp(k)=pars(k)+myeps;
  disimkerntemp=kernExpandParam(disimkern,[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  simkerntemp=kernExpandParam(simkern,[parstemp(1) parstemp(2) parstemp(3)]);
  Ktemp2=disimXsimKernCompute(disimkerntemp,simkerntemp,model.t);

  diffkerngradient2b(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));
end;
modelkerngradient2b
diffkerngradient2b

pause

end;



if 0,

covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
%covGrad=ones(length(model.t)*2);
modelkerngradient3=kernGradient(model.kern,model.t,covGrad);
[pars,nams]=kernExtractParam(model.kern);

for k=1:length(pars),
  parstemp=pars;
  myeps=1e-5;
  
  parstemp(k)=pars(k)-myeps;
  kerntemp=kernExpandParam(model.kern,parstemp);
  Ktemp1=kernCompute(kerntemp,model.t);
    
  parstemp(k)=pars(k)+myeps;
  kerntemp=kernExpandParam(model.kern,parstemp);
  Ktemp2=kernCompute(kerntemp,model.t);

  kernderiv=(Ktemp2-Ktemp1)/(2*myeps);
  diffkerngradient3(k)=sum(sum(kernderiv.*covGrad));
end;
modelkerngradient3
diffkerngradient3

pause





kern=model.kern;
pars=[kern.comp{1}.comp{1}.decay kern.comp{1}.comp{1}.inverseWidth kern.comp{1}.comp{1}.variance kern.comp{1}.comp{2}.decay kern.comp{1}.comp{2}.variance kern.comp{1}.comp{2}.rbf_variance kern.comp{2}.comp{1}.variance kern.comp{2}.comp{2}.variance];
%covGrad=ones(length(model.t)*2);
covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;

modelkerngradient4=cmpndKernGradient(kern,model.t,covGrad);
kfactors=kernFactors(kern,'gradfact');

for k=1:length(pars),
  parstemp=pars;
  myeps=1e-5;

  fprintf(1,'Here A\n');
  parstemp(k)=pars(k)-myeps;
  tempkern=kern;
  tempkern.comp{1}.comp{1}=simKernExpandParam(kern.comp{1}.comp{1},[parstemp(1) parstemp(2) parstemp(3)]);
  tempkern.comp{1}.comp{2}=disimKernExpandParam(kern.comp{1}.comp{2},[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  tempkern.comp{2}.comp{1}=whiteKernExpandParam(kern.comp{2}.comp{1},[parstemp(7)]);
  tempkern.comp{2}.comp{2}=whiteKernExpandParam(kern.comp{2}.comp{2},[parstemp(8)]);
%  tempkern.comp{1}.comp{1}
%  tempkern.comp{1}.comp{2}
%  pause
  Ktemp1=kernCompute(tempkern,model.t);
    
  fprintf(1,'Here B\n');
  parstemp(k)=pars(k)+myeps;
  tempkern=kern;
  tempkern.comp{1}.comp{1}=simKernExpandParam(kern.comp{1}.comp{1},[parstemp(1) parstemp(2) parstemp(3)]);
  tempkern.comp{1}.comp{2}=disimKernExpandParam(kern.comp{1}.comp{2},[parstemp(1) parstemp(2) parstemp(3) parstemp(4) parstemp(5) parstemp(6)]);
  tempkern.comp{2}.comp{1}=whiteKernExpandParam(kern.comp{2}.comp{1},[parstemp(7)]);
  tempkern.comp{2}.comp{2}=whiteKernExpandParam(kern.comp{2}.comp{2},[parstemp(8)]);
  Ktemp2=kernCompute(tempkern,model.t);

  diffkerngradient4(k)=sum(sum((Ktemp2-Ktemp1)/(2*myeps).*covGrad));

end;

factors(1)=sigmoidabTransform(pars(1),'gradfact',kern.comp{1}.comp{1}.transforms(1).transformsettings);
factors(2)=sigmoidabTransform(pars(2),'gradfact',kern.comp{1}.comp{1}.transforms(2).transformsettings);
factors(3)=sigmoidabTransform(pars(3),'gradfact',kern.comp{1}.comp{1}.transforms(3).transformsettings);
factors(4)=sigmoidabTransform(pars(4),'gradfact',kern.comp{1}.comp{2}.transforms(4).transformsettings);
factors(5)=sigmoidabTransform(pars(5),'gradfact',kern.comp{1}.comp{2}.transforms(5).transformsettings);
factors(6)=sigmoidabTransform(pars(1),'gradfact',kern.comp{1}.comp{2}.transforms(6).transformsettings);
factors(7)=sigmoidabTransform(pars(1),'gradfact',kern.comp{2}.comp{1}.transforms(1).transformsettings);
factors(8)=sigmoidabTransform(pars(1),'gradfact',kern.comp{2}.comp{2}.transforms(1).transformsettings);
diffkerngradient4=diffkerngradient4.*factors;

modelkerngradient4
diffkerngradient4

pause


end;



if 0,

  % Compute the part of the gradient where the mean is held fixed
  
  covGrad = -model.invK + model.invK*model.m*model.m'*model.invK;
  covGrad = 0.5*covGrad;
  kern=model.kern;  
  modelkerngradient5=kernGradient(kern,model.t,covGrad);
  modelkerngradient5b=gpdisimLogLikeGradients(model,1,0);
  
  [pars,nams]=gpdisimExtractParam(model);
  diffkerngradient=zeros(size(pars));
  myeps=1e-5;
  for k=1:length(pars),
    parstemp=pars;
    
    parstemp(k)=pars(k)-myeps;
    modeltemp=gpdisimExpandParam(model,parstemp,1,0);
    valtemp1=gpdisimLogLikelihood(modeltemp);
    
    parstemp(k)=pars(k)+myeps;
    modeltemp=gpdisimExpandParam(model,parstemp,1,0);
    valtemp2=gpdisimLogLikelihood(modeltemp);
    
    diffkerngradient5(k)=(valtemp2-valtemp1)/(2*myeps);
  end;
  
  modelkerngradient5
  modelkerngradient5b
  diffkerngradient5
  
end;
  

if 0,
  [pars,nams]=gpdisimExtractParam(model);
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
  end;
  
  modeltotalgradient
  difftotalgradient
end;



% df/da = (df/dy) * (dy/dx)
% y=A+(B-A)/(1+exp(-x))   
% ----> dy/dx 
% = -(B-A)/(1+exp(-x))^2*(-1)*exp(-x) 
% = (B-A)*(1/(1+exp(-x)))*(exp(-x)/(1+exp(-x)))
% = (B-A)*((y-A)/(B-A))*(1-(y-A)/(B-A))
% = (y-A)*(B-A - (y-A))/(B-A)
% = (y-A)*(B-y)/(B-A)
