function [loglikelihood,loglikelihoodsperpoint,ll1perpoint,ll2perpoint]=computeLeaveOneOutLikelihood(model);

[pars,nams]=gpdisimExtractParam(model);
nt=length(model.t);
vals1=model.y(1:nt);
vals2=model.y(nt+1:2*nt);

inversewidth_index=2;
inversewidth=model.kern.comp{1}.comp{1}.inverseWidth;
lengthscale=sqrt(1/(2*inversewidth));



loglikelihoodsperpoint=[];
for k=1:length(model.t),
  tempvals1=vals1([(1:k-1) (k+1:nt)]);
  tempvals2=vals2([(1:k-1) (k+1:nt)]);
  temptimes=model.t([(1:k-1) (k+1:nt)]);
  tempmodel=createSimDisim(temptimes,tempvals1,tempvals2,lengthscale,1);
  %tempmodel.fix
  %pause
  
  tempmodel=gpdisimExpandParam(tempmodel,pars);
  [pars2,nams2]=gpdisimExtractParam(model);
%  pars2
%  nams2
  
  [temppriormeans,tempmeans,covmatrix]=gpasimTemp4Predict(tempmodel,model.t(k),1);
%  covmatrix
%  pause

  
  dim=2;
  pred_loglikelihood=-(dim/2)*log(2*pi) -(1/2)*log(det(covmatrix)) ...
      - (1/2)*([vals1(k);vals2(k)]-tempmeans)'*inv(covmatrix)*([vals1(k);vals2(k)]-tempmeans);
  loglikelihoodsperpoint(k)=pred_loglikelihood;

  dim=1;
  pred_loglikelihood=-(dim/2)*log(2*pi) -(1/2)*log(det(covmatrix(1,1))) ...
      - (1/2)*([vals1(k)]-tempmeans(1))'*inv(covmatrix(1,1))*([vals1(k)]-tempmeans(1));
  ll1perpoint(k)=pred_loglikelihood;

  dim=1;
  pred_loglikelihood=-(dim/2)*log(2*pi) -(1/2)*log(det(covmatrix(2,2))) ...
      - (1/2)*([vals2(k)]-tempmeans(2))'*inv(covmatrix(2,2))*([vals2(k)]-tempmeans(2));
  ll2perpoint(k)=pred_loglikelihood;
end;
%ll1perpoint
%ll2perpoint
%pause
loglikelihood=sum(loglikelihoodsperpoint);
