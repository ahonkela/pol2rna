function geneprofiles = make_geneprofiles_v2(utrlocations,tempbed);

%
% Find BEDGRAPH lines corresponding to each of the known 24
% chromosomes (1-22, X, Y)
%
chromosomeI=cell(24,1);
for k=1:24,
  chromosomeI{k}=find(tempbed(:,1)==k);
end;

binsize=200;
ngenes=size(utrlocations,1);
geneprofiles=cell(ngenes,1);

%
% Process each gene in turn, and compute its groseq profile.
% The code assumes that areas described in different lines in the BEDGRAPH file are nonoverlapping.

for k=1:ngenes,
    if mod(k,100)==0,
        k
    end;
    
    % get gene information (chromosome index, gene area). Note that
    % UTR3 area is not used for computing profiles, profile is computed
    % over the whole gene!
    chro=utrlocations(k,7);
    genestart=utrlocations(k,5);
    geneend=utrlocations(k,6);
    strand=utrlocations(k,4);
    
    % length of the gene in basepairs
    genelength=geneend-genestart+1;
    
    % length of the gene in 200bp bins
    binlength=ceil(genelength/binsize);

    % We first compute the profile in detail over each basepair,
    % and later scale it down to a profile over 200bp bins
    geneprofile=zeros(genelength,1);

    % find BED lines that are in the gene area
    Ichro=chromosomeI{chro};
    I=find((tempbed(Ichro,3)>=genestart)&(tempbed(Ichro,2)<=geneend));
    if length(I)>0,
        %truncate regions to gene area
        I=Ichro(I);
        templines=tempbed(I,:);    
        I2=find(templines(:,2)<genestart);templines(I2,2)=genestart;
        I2=find(templines(:,3)>geneend);templines(I2,3)=geneend;
        %fill regions into gene profile
        templines(:,2)=templines(:,2)-genestart+1;
        templines(:,3)=templines(:,3)-genestart+1;
        for m=1:length(I),
          geneprofile(templines(m,2):templines(m,3))=templines(m,4);
        end;
    end;
   
    % scale down profile to 200bp bins starting from transcription
    % end. UTR3 region will be at the few first elements of the bin
    % profile!

    % For genes on the positive strand, the profile needs to be
    % reversed so that the UTR3 area is at the start. For genes on
    % the negative strand, this is already the case.
    if strand==1,
      geneprofile=geneprofile(end:-1:1);
    end;

    % scale down the profile to a profile over the 200bp bins.
    % Each bin reports the mean activity over basepairs of that
    % bin. The last bin is handled separately since it may be
    % shorter than 200bp.
    
    binprofile=zeros(binlength,1);
    for m=1:binlength-1,
      binprofile(m)=mean(geneprofile((m-1)*binsize + [1:binsize]));
    end;
    binprofile(binlength)=mean(geneprofile((binlength-1)*binsize+1:end));

    geneprofiles{k}=binprofile;    
end;
