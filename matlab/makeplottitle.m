function plottitle=makeplottitle(model,loglik,ensemblid);

ndigits=2;
Bstr=smallprint(model.B(1),ndigits);
Dstr=smallprint(model.D(1),ndigits);
Sstr=smallprint(model.S(1),ndigits);
delaystr=smallprint(model.delay(1),ndigits);
MuPolstr=smallprint(model.simMean,ndigits);
Rna0str=smallprint(model.disimStartMean,ndigits);
nvarPolstr=smallprint(model.kern.comp{2}.comp{1}.variance,ndigits);

plottitle=sprintf('ENSG%011d, LL=%f, B=%s D=%s\nS=%s delay=%s mu_{Pol2}=%s Rna_0=%s nvar_{Pol2}=%s\n',...
  ensemblid,... 
  loglik, ...
  Bstr, Dstr, Sstr, delaystr, MuPolstr, Rna0str, nvarPolstr);
