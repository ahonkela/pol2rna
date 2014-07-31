function [bestmodel,joint_ll]=findbestDisimModel(timevector,dataVals1,dataVals2,lengthscale,ninits,maxiters);



% joint model
%for k=2:6,
for k=0:0,
  % leave one time point out of POL2 and RNA
  temptimes=timevector([(1:k-1) (k+1:length(timevector))]);
  tempvals1=dataVals1([(1:k-1) (k+1:length(dataVals1))]);
  tempvals2=dataVals2([(1:k-1) (k+1:length(dataVals2))]);

  % normalization of POL2
  tempvals1=tempvals1-dataVals1(1);
  tempvals1=tempvals1/sqrt(var(dataVals1));

  if ~isempty(tempvals2),
    % normalization of RNA
    tempvals2=tempvals2-dataVals2(1);
    tempvals2=tempvals2/sqrt(var(dataVals2));
  end;

  joint_ll=-inf;
  for initializationtype=1:10,
    jointmodel=createSimDisim(temptimes,tempvals1,tempvals2,lengthscale,initializationtype);

    [pars,nams]=gpdisimExtractParam(jointmodel);
    %ninits=5;
    for m=1:ninits,
      if m==1,
        temppars=pars;
      else    
        temppars=pars+7*randn(size(pars));
      end;
      tempmodel=gpdisimExpandParam(jointmodel,temppars);
      try
        tempmodelb = gpdisimOptimise(tempmodel,1,maxiters);
      catch
        tempmodelb=tempmodel;
      end;
    
      % templl=gpdisimLogLikelihood(tempmodelb);
      [loglikelihood,loglikelihoodsperpoint]=computeLeaveOneOutLikelihood(tempmodelb);
      templl=loglikelihood;
      
      if templl>joint_ll,
        jointmodelb=tempmodelb;
        joint_ll=templl;
      end;  
    end;
  end;
  
  %jointmodelb=gpdisimOptimise(jointmodel,1,maxiters);
  %joint_ll=gpdisimLogLikelihood(jointmodelb);

  if k>0,
    % predict values at the missing time point
    [temppriormeans,tempmeans,covmatrix]=gpasimTemp4Predict(gpsim3model,timevector(k));
    preddim=2;
    tempm=[(dataVals1(k)-mean(dataVals1))/sqrt(var(dataVals1));...
           (dataVals2(k)-mean(dataVals2))/sqrt(var(dataVals2))]-tempmeans;
    pred_ll= -preddim*log(2*pi) - log(det(covmatrix)) - tempm'*inv(covmatrix)*tempm;
  end;
end;

bestmodel=jointmodelb;
