function [mypoldata,myrnadata,mytimes]=generateDataFromModel(model,times);

model=model;
[pars,nams]=gpdisimExtractParam(model);

model.t=times;
model.yvar=0*[times;times];
model.y=model.yvar;
model=gpdisimExpandParam(model,pars);



mu=model.mu;
K=model.K;


tempdata=mvnrnd(mu,K)';
%size(tempdata)
mypoldata=tempdata(1:length(times));
myrnadata=tempdata(length(times)+1:2*length(times));
mytimes=times;