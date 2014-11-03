---------------------------------------------------------------------
% Note that the codes presented here are given only for reference.
% If you intend to use these codes, please first consider contacting
% hande.topa@aalto.fi for more detailed instructions.
---------------------------------------------------------------------


1. Match the transcripts with genes:

transcript_file='Homo_sapiens.GRCh37.69.all.tr';
indices=getIndices(transcript_file); 
u_indices=unique(indices);

2. Collect the BitSeq file names in a cell:

time_points={'t0000','t0005','t0010','t0020','t0040','t0080','t0160','t0320','t0640','t1280'};


3. Create mu (for mean log rpkm gene expression levels) and V (for technical variances of log-rpkm gene expression levels) matrices:

mu=[];
V=[];

4. Read the data files and get the required matrices for biological variance estimation and GP analysis:

for t=1:10
	filename=[time_point{t},'.rpkm'];
	f0=fopen(filename)
	t0=textscan(f0,repmat('%f',1,500)'CommentStyle','#');
	fclose(f0)
	tr_mcmc_samples=cell2mat(t0);

	gene_mcmc_samples=zeros(length(u_indices),500);

	for i=1:length(u_indices)
		f=find(indices==i);
		gene_mcmc_samples(i,:)=sum(tr_mcmc_samples(f,:),1);
	end

	if t==1
		t1_mcmc=gene_mcmc_samples;
	else if t==2
		t2_mcmc=gene_mcmc_samples;
	else if t==3
		t3_mcmc=_gene_mcmc_samples;
	end
	
	mu=[mu mean(log(gene_mcmc_samples),2)];
	V=[V var(log(gene_mcmc_samples),0,2)];

end

5. Estimate the biological variances for gene groups:

[m_groups biol_var]=estimate_biolvar(t1_mcmc,t2_mcmc,t3_mcmc);

6. Run GP analysis to get the Bayes factors, which will be saved in the first column of the output file named 'GP_summary':

GPanalysis(mu,V,m_groups,biol_var);


