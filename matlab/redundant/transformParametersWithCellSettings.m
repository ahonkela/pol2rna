function transformedparams=transformParametersWithSettings(params,transforminfo,transformtype);

if isempty(transformtype),
  transformtype='atox';
end;


transformedparams=params;
for k=1:length(transforminfo),
    tempindex=k;
    temprange=transforminfo{k};
    if strcmp(transformtype,'atox')==1,
      transformedparams(tempindex)=sigmoidabTransform(params(tempindex),'atox',temprange); 
    elseif strcmp(transformtype,'xtoa')==1,
      transformedparams(tempindex)=sigmoidabTransform(params(tempindex),'xtoa',temprange);       
    elseif strcmp(transformtype,'gradfact')==1,
      transformedparams(tempindex)=sigmoidabTransform(params(tempindex),'gradfact',temprange);       
    end;     
end;
