% Create comparison plots for comparing the new RNA and PolII 
% data (March 2012) to the old data in Tigrebrowser at 
% http://users.ics.tkk.fi/jtpelto/tigrebrowser/tigreBrowser.cgi

cd /share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed
load pol2_for_matti_ver3.mat  % Needed for 'bininfo' structure
load pol2_summaryseries_2012_03.mat % provides pol2_summaryseries

addpath /share/mi/workspace/jtpelto/synergy/mlprojects/mlprojects/pol2rnaseq/matlab
rnafile = read_stringfile('/share/synergy/data_summaries/new_rna_counts.txt',[],[' ' 10 13 0]);
rna_summaryseries=nan*ones(size(bininfo,1),10);
for k=2:length(rnafile),
  if mod(k,1000) == 0,
    k
  end;
  ensgname = rnafile{k}{1};
  ensgid = str2double(ensgname(5:end));
  I = find(bininfo(:,5)==ensgid);
  if length(I) > 0,
    for l=1:10,
      rna_summaryseries(I,l) = str2double(rnafile{k}{1+l});
    end;
  end;
end;


% normalization of POL2
pol2_summaryseries_normalized = pol2_summaryseries;
for k=1:size(bininfo,1),
  dataVals1 = pol2_summaryseries(k,:);
  dataVals1=dataVals1-min(dataVals1);
  if mean(dataVals1.^2)>0,
    dataVals1=dataVals1/sqrt(mean(dataVals1.^2));
  end;
  pol2_summaryseries_normalized(k,:) = dataVals1;  
end;


% RNA normalization by overall read count per time point
rna_summaryseries_normalized = rna_summaryseries;
I = find(isnan(rna_summaryseries(:,1))==0);
for k=1:10,
  overallreadcount = sum(rna_summaryseries(I,k));
  rna_summaryseries_normalized(:,k) = rna_summaryseries(:,k)/overallreadcount;
end;

% normalization of RNA per time series
for k=1:size(bininfo,1),
  dataVals2 = rna_summaryseries_normalized(k,:);
  if mean(dataVals2.^2)>0,
    normalizationfactor=sqrt(mean(dataVals2.^2));
    dataVals2=dataVals2/normalizationfactor;
  end;
  rna_summaryseries_normalized(k,:) = dataVals2;
end;

% tigrebrowser plot names
tigrebrowserplots = read_stringfile('/share/mi/workspace/jtpelto/synergy/synergy_data/PolII/processed/polrnaconnectionresults.txt',[],[' ' 10 13 0]);

for k=2:length(tigrebrowserplots),
  browserindices(k-1) = str2double(tigrebrowserplots{k}{1}(5:end));
end;

% strongpol2 = find(sum(pol2_summaryseries,2) > 10*(5e5));
% for k=1:length(strongpol2),
%   isinbrowser(k) = sum(bininfo(strongpol2(k),5)==browserindices);
% end;
% I = find(isinbrowser==1);
% isinbrowser_all = zeros(size(bininfo,1),1);
% isinbrowser_all(strongpol2(I)) = 1;

isinbrowser_all = zeros(size(bininfo,1),1);
for k=1:size(pol2_summaryseries,1),
  if (sum(pol2_summaryseries_normalized(k,:))>0) ...
      &&(sum(rna_summaryseries_normalized(k,:))>0) ...
      &&(sum(browserindices==bininfo(k,5))>0),
    isinbrowser_all(k) = 1;
  end;
end;

f = fopen('newdataresults.txt','w');
I = find(isinbrowser_all == 1);
fprintf(f,'ENSID Pol2Avg RNAAvg\n');
for k=1:length(I),
  geneindex=I(k);
  genelength=bininfo(geneindex,3)-bininfo(geneindex,2)+1;
  pol2avg = sum(pol2_summaryseries(geneindex,:))/(10);
  rnaavg = sum(rna_summaryseries(geneindex,:))/(10*genelength);
  fprintf(f, 'ENSG%011d %e %e\n', bininfo(geneindex,5), pol2avg, rnaavg);
  clf;
  subplot(2,1,1);
  h=plot(sqrt([0 5 10 20 40 80 160 320 640 1280]),pol2_summaryseries_normalized(geneindex,:),'r-');
  set(h,'LineWidth',2);
  h=gca;
  set(h,'XTick',sqrt([0 5 10 20 40 80 160 320 640 1280]));
  set(h,'XTickLabel',{'0','5','10','20','40','80','160','320','640','1280'});  
  set(h,'XLim',[0 sqrt(1280)]);
  grid on;
  subplot(2,1,2);
  h = plot(sqrt([0 5 10 20 40 80 160 320 640 1280]),rna_summaryseries_normalized(geneindex,:),'g-');
  set(h,'LineWidth',2);
  h=gca;
  set(h,'XTick',sqrt([0 5 10 20 40 80 160 320 640 1280]));
  set(h,'XTickLabel',{'0','5','10','20','40','80','160','320','640','1280'});  
  set(h,'XLim',[0 sqrt(1280)]);
  grid on;
  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7 6])
  print('-dpng','-r100',sprintf('ENSG%011d_pol2rna_new.png',bininfo(geneindex,5)));
end;
fclose(f);

