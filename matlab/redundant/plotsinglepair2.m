function plotsinglepair2(polseries,rnaseries,geneindices,geneinfo,k);

  serieslength=min([size(polseries,2) size(rnaseries,2)]);

  timepoints=[0 5 10 20 40 80 160 320 640 1280];
  timepoints=timepoints(1:serieslength);

  pol=polseries(k,1:serieslength);
  rna=rnaseries(k,1:serieslength);

  geneindex=geneindices(k);
  gene_startpoint = geneinfo{geneindex}{3}(2:end-1);
  gene_endpoint = geneinfo{geneindex}{4}(2:end-1);
  genename = geneinfo{geneindex}{1}(2:end-1);
  genechr = geneinfo{geneindex}{2}(2:end-1);


  clf;
  subplot(2,1,1);
  h=plot(timepoints,pol,'r-');
  axis tight;
  set(h,'LineWidth',2);

  pol2title = sprintf('POL2 over the gene, ');
  title(pol2title);
  
  subplot(2,1,2);
  h=plot(timepoints,rna,'b-');
  axis tight;
  set(h,'LineWidth',2);


 % tempcorr=corr([pol2peakheights' rna']);  
 % pearsoncorr=tempcorr(1,2);
if (var(pol)>0) && (var(rna)>0), 
pearsoncorr=((pol-mean(pol))*(rna-mean(rna))')/serieslength/sqrt(var(pol,1)*var(rna,1));
else 
pearsoncorr=0;
end;

 [tempcorr,pvals]=corr([pol' rna']);
 pearsoncorr2=tempcorr(1,2);
 pval=pvals(1,2);


% compute average squared deviation of RNA from zero
  rnanonzero=find(rna~=0);
  if length(rnanonzero) > 0,
    % assume that all exactly zero scores are missing values, compensate
    rnavar = sum(rna(rnanonzero).^2)*serieslength/length(rnanonzero);
  else
    rnavar=1;
  end;
   
  % compute average squared deviation of POL2 from zero
  pol2nonzero=find(pol~=0);
  if length(pol2nonzero) > 0,
    % assume that all exactly zero scores are missing values, compensate
    pol2var = sum(pol(pol2nonzero).^2)*serieslength/length(pol2nonzero);
  else
    pol2var = 1;
  end;

  % compute unnormalized correlation (just sum of products)
  sum_rna_to_pol2 = pol*rna';

  % compute the "assumed zero-mean, missing value compensated" correlation
  zeromeancorrelation = sum_rna_to_pol2/sqrt(pol2var*rnavar);


  rnatitle=sprintf('Gene %s (%s: %s - %s)\nRNA based expression\nPearson corr. %f (pval %f)\nzero-mean corr. %f', ...
    genename, genechr, gene_startpoint, gene_endpoint,pearsoncorr2,pval,zeromeancorrelation);
  title(rnatitle);
  

