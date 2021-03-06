function [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createGeneGPModels(timevector,dataVals1,dataVals2,lengthscale,maxiters,ninits);


% Naive log-likelihood of POL2
if mean(dataVals1.^2)>0,
  datatemp=dataVals1/sqrt(mean(dataVals1.^2));
  dim = size(datatemp, 1);
  tempK=eye(dim)*var(datatemp);
  tempm=datatemp-mean(datatemp);
  ll1 = -(dim/2)*log(2*pi) - (1/2)*log(det(tempK)) - (1/2)*tempm'*inv(tempK)*tempm;
else
  ll1=inf;
end;

% Naive log-likelihood of RNA
if mean(dataVals2.^2)>0,
  datatemp=dataVals2/sqrt(mean(dataVals2.^2));
  dim = size(datatemp, 1);
  tempK=eye(dim)*var(datatemp);
  tempm=datatemp-mean(datatemp);
  ll2 = -(dim/2)*log(2*pi) - (1/2)*log(det(tempK)) - (1/2)*tempm'*inv(tempK)*tempm;
else
  ll2=inf;
end;

% Naive log-likelihood of independently modeled POL2 and RNA
naive_ll=ll1+ll2
%pause

create_rbfmodels=0;

if create_rbfmodels==1,
  %--------------
  % Independent ARBF models for both POL2 and RNA
  %--------------
  pol2_ll=-inf;
  for initializationtype=1:ninits,
    pol2model = createRBFModel(timevector,dataVals1,lengthscale,initializationtype);
    [pars,nams]=gpdisimExtractParam(pol2model);
    %ninits=5;
    temppars=pars;
    temppars
    tempmodel=gpdisimExpandParam(pol2model,temppars);
    fprintf(1,'POL2 rbf model numGenes %d\n',tempmodel.numGenes);
    tempmodelb = gpdisimOptimise(tempmodel,1,maxiters);
    templl=gpdisimLogLikelihood(tempmodelb)
    % pause
    if templl>pol2_ll,
      pol2modelb=tempmodelb;
      pol2_ll=templl;
    end;  
  end;
  
  
  rna_ll=-inf;
  for initializationtype=1:ninits,
    rnamodel = createRBFModel(timevector,dataVals2,lengthscale,initializationtype);
    [pars,nams]=gpdisimExtractParam(rnamodel);
    %ninits=5;
    temppars=pars;
    tempmodel=gpdisimExpandParam(rnamodel,temppars);
    tempmodelb = gpdisimOptimise(tempmodel,1,maxiters);
    templl=gpdisimLogLikelihood(tempmodelb);
    if templl>rna_ll,
      rnamodelb=tempmodelb;
      rna_ll=templl;
    end;  
  end;

  rbf_ll=pol2_ll+rna_ll
  %pause
else
  pol2modelb=nan;
  rnamodelb=nan;
  rbf_ll=nan;
end;




if 0,
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
    [jointmodel,temptransforminfo]=createSimDisim(temptimes,tempvals1,tempvals2,lengthscale,initializationtype);
    [pars,nams]=gpdisimExtractParam(jointmodel);
    %ninits=5;
    %for m=1:ninits,
    %  if m==1,
    %    temppars=pars;
    %  else    
    %    temppars=pars+7*randn(size(pars));
    %  end;
    %  tempmodel=gpdisimExpandParam(jointmodel,temppars);
    %  try
    %    tempmodelb = gpdisimOptimise(tempmodel,1,maxiters);
    %  catch
    %    tempmodelb=tempmodel;
    %  end;

    try
      tempmodelb = gpdisimOptimise(jointmodel,1,maxiters);
    catch
      tempmodelb=jointmodel;
    end;

    
    
    % templl=gpdisimLogLikelihood(tempmodelb);
    [loglikelihood,loglikelihoodsperpoint]=computeLeaveOneOutLikelihood(tempmodelb);
    templl=loglikelihood
    pause      
    if templl>joint_ll,
      jointmodelb=tempmodelb;
      jointtransforminfo=temptransforminfo;
      joint_ll=templl;
    end;  
    %end;
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
end;




joint_ll=-inf;
for initializationtype=1:ninits,
  try
    fprintf(1,'Trying to create initialization %d for SIM-DISIM\n',initializationtype);
    [jointmodel,temptransforminfo]=createSimDisim(timevector,dataVals1,dataVals2,lengthscale,initializationtype);
    fprintf(1,'Trying to create initialization %d for SIM-DISIM, step2\n',initializationtype);
    init_ll=gpdisimLogLikelihood(jointmodel);
    %  init_ll
    %pause
  
    fprintf(1,'Trying to create initialization %d for SIM-DISIM, step3\n',initializationtype);
    [pars,nams]=gpdisimExtractParam(jointmodel)
    % jointmodel.disimStartMean
    % jointmodel.y
    % jointmodel.mu
    % jointmodel.m
    
    %plotpredictions(jointmodel,[0:5:1280]',2,1,1,'exampletitle');
    %drawnow;
    %pause

    fprintf(1,'Trying to optimize from initialization %d for SIM-DISIM\n',initializationtype);
    try
      tempmodelb = gpdisimOptimise(jointmodel,1,maxiters);
    catch
      tempmodelb=jointmodel;
    end;
        
    fprintf(1,'Trying to evaluate model trained from initialization %d for SIM-DISIM\n',initializationtype);
    templl=gpdisimLogLikelihood(tempmodelb);
    % [templl,loglikelihoodsperpoint]=computeLeaveOneOutLikelihood(tempmodelb);
    
    %plotpredictions(tempmodelb,[0:5:1280]',2,1,1,'exampletitle');
    fprintf(1,'Initialization %d: init_ll %f, optimized %f\n', initializationtype, init_ll, templl);
    %pause      
    fprintf(1,'Done processing initialization %d for SIM-DISIM\n',initializationtype);
  catch
    templl=-inf;
  end;
    
  if templl>joint_ll,
    jointmodelb=tempmodelb;
    jointtransforminfo=temptransforminfo;
    joint_ll=templl;
  end;  
end;




% onko SIM-DISIM koodissa gradienttiongelma? Ei vastaa empiirista
% gradienttia taysin!




