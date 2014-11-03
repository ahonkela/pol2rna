function[I]=getIndices(transcript_file)

f0=fopen(transcript_file)
t0=textscan(f0,'%s%s%f%f','CommentStyle','#')
fclose(f0)
all_genes=t0{1};
trs=t0{2};


genes=unique(all_genes,'stable');
N=length(genes); % number of genes
M=length(trs); % number of transcripts

I=zeros(M,1);

for i=1:N
	gene_name=genes(i);
	tr_inds=find(strcmp(gene_name,all_genes)==1);
        I(tr_inds)=i;
end

end



