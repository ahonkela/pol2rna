function [difftotalgradient]=testKernelTotalGradient(kernel,t,covGrad);

[pars,nams]=kernExtractParam(kernel);
kernel=kernExpandParam(kernel,pars);


[simdiffgradient,simgradient]=testKernelGradient(kernel.comp{1}.comp{1},t,covGrad(1:9,1:9));

[disimdiffgradient,disimgradient]=testKernelGradient(kernel.comp{1}.comp{2},t,covGrad(10:18,10:18));

[diffcrossgradient,crossgradient]=testDisimSimCrossGradient(kernel.comp{1}.comp{2},kernel.comp{1}.comp{1},t,covGrad(10:18,1:9));

simdiffgradient
%  decay       invwidth      variance
%  2.4812e+06  -5.4911e+04   5.8549e+03
%   2.4662e+06  -5.4915e+04   5.8549e+03  analytic


disimdiffgradient
%  simdecay    invwidth      simvariance  disimdecay  disimvariance rbfvariance
%  2.3518e+06  -5.4857e+04   5.8509e+03  -1.4209e+04   5.8509e+03   2.9254e+03 
%  2.3554e+06  -5.4850e+04   5.8509e+03  -1.4167e+04   5.8509e+03   2.9254e+03  analytic

diffcrossgradient
%  -2.4602e+06   5.4971e+04  -5.8553e+03   7.0975e+03  -2.9277e+03  -2.9277e+03
%  -2.5035e+06    5.4960e+04  -5.8553e+03   7.0856e+03  -2.9277e+03  -2.9277e+03 analytic




totalgradient=diffcrossgradient*2+disimdiffgradient+[simdiffgradient 0 0 0]

% -8.7463e+04   1.7376e+02  -4.8954e+00  -1.4341e+01  -4.4315e+00  -2.9299e+03



%analytic compgradient
%  -1.8525e+05   1.5554e+02  -4.8952e+00   4.5672e+00  -4.4314e+00  -2.9299e+03

% diff testKerneCompGradient
% -8.7463e+04   1.7376e+02  -4.8954e+00  -1.4341e+01  -4.4315e+00  -2.9299e+03   0.0000e+00   0.0000e+00

%multiKernGradient for comp1 before paramgroups multiplication
   simdecay    invwidth      simvar      simdecay      invwidth     simvar       disimdecay   disimvar    rbfvar
   2.4662e+06  -5.4915e+04   5.8549e+03  -2.6515e+06   5.5071e+04  -5.8597e+03   4.5672e+00  -4.4314e+00  -2.9299e+03

