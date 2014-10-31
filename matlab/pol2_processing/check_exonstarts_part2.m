% Compare each gene start to the area in bininfo
nproblemgenes=0;
gene_originallengths=nan*ones(size(geneids));
gene_strands=nan*ones(size(geneids));
n_missinggenes=0;
for k=1:length(geneids),
  if isnan(genestarts(k))==0,
    Itemp=find(bininfo(:,5)==geneids(k));
    if length(Itemp==1),
      % Gene is found in bininfo structure
      if bininfo(Itemp,6)==1,  
        % gene orientation is the forward strand
        length1 = bininfo(Itemp,3)-bininfo(Itemp,2)+1; % length from gene start to end
        length2 = bininfo(Itemp,3)-genestarts(k)+1;  % length from active exon to gene end
        extendedstart = genestarts(k);
        extendedend = bininfo(Itemp,3);
      else
        % gene orientation is the reverse strand
        length1 = bininfo(Itemp,3)-bininfo(Itemp,2)+1; % length from gene start to end
        length2 = -(bininfo(Itemp,2)-geneends(k))+1;  % length from active exon to gene end
        extendedstart = bininfo(Itemp,2);
        extendedend = geneends(k);
      end;
      gene_originallengths(k)=length1;
      gene_strands(k)=bininfo(Itemp,6);
      if length2<=0.5*length1,
        nproblemgenes=nproblemgenes+1;
        problemgenes(nproblemgenes)=geneids(k);
        problemgene_newstarts(nproblemgenes)=genestarts(k);
        problemgene_newends(nproblemgenes)=geneends(k);
        problemgene_newlengths(nproblemgenes)=length2;
        problemgene_strands(nproblemgenes)=bininfo(Itemp,6);
        problemgene_origlengths(nproblemgenes)=length1;
        problemgene_extendedstarts(nproblemgenes)=extendedstart;
        problemgene_extendedends(nproblemgenes)=extendedend;
      end;
    else
      % Gene is found in bininfo structure
      n_missinggenes=n_missinggenes+1;
      missinggenes(n_missinggenes)=geneids(k);
    end;
  end;
end;



