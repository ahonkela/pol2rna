function [jointmodelb,jointtransforminfo,pol2modelb,rnamodelb,naive_ll,rbf_ll,joint_ll]=createNdGeneGPModels_longdelay(timevector,dataVals1,dataVals2,lengthscale,maxiters,ninits,parameterranges,use_fixedrnavariance);


% Naive log-likelihood of POL2
datatemp=dataVals1-min(dataVals1);  
if mean(datatemp.^2)>0,
  datatemp=datatemp/sqrt(mean(datatemp.^2));
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







joint_ll=-inf;
for initializationtype=1:ninits,
  try
    fprintf(1,'Trying to create initialization %d for SIM-DISIM\n',initializationtype);
    [jointmodel,temptransforminfo]=createNdSimDisim_longdelay(timevector,dataVals1,dataVals2,lengthscale,initializationtype,parameterranges,use_fixedrnavariance);
    fprintf(1,'Trying to create initialization %d for SIM-DISIM, step2\n',initializationtype);
    init_ll=gpnddisimLogLikelihood(jointmodel);
    %  init_ll
    %pause
  
    fprintf(1,'Trying to create initialization %d for SIM-DISIM, step3\n',initializationtype);
    [pars,nams]=gpnddisimExtractParam(jointmodel)
    % jointmodel.disimStartMean
    % jointmodel.y
    % jointmodel.mu
    % jointmodel.m
    
    %plotpredictions(jointmodel,[0:5:1280]',2,1,1,'exampletitle');
    %drawnow;
    %pause

    %pars
    %nams
    %pause
    %testKernelGradient(jointmodel.kern,jointmodel.t,jointmodel.invK)
    %pause
    
    fprintf(1,'Trying to optimize from initialization %d for SIM-DISIM\n',initializationtype);
    try
      tempmodelb = gpnddisimOptimise(jointmodel,1,maxiters);
    catch
      tempmodelb=jointmodel;
    end;
        
    fprintf(1,'Trying to evaluate model trained from initialization %d for SIM-DISIM\n',initializationtype);
    templl=gpnddisimLogLikelihood(tempmodelb);
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




