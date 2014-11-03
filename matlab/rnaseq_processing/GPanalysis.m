% This code should be run under the same directory where gptk matlab functions are placed.
% Inputs:
% - mu: 39071-by-10 matrix which contains mean log gene expression levels (log-rpkm)
% - V: 39071-by-10 matrix which contains the technical variances of log gene expression levels

function[]=GPanalysis(mu,V,m_groups,biol_var)

t=[0 5 10 20 40 80 160 320 640 1280]';
x=sqrt(t);

N=length(mu); % number of genes: 39071
BF=zeros(N,1); % Bayes factors
likli0=zeros(N,1); % likelihood of null model
likli1=zeros(N,1); % likelihood of alternative model
LS=zeros(N,1); % length-scale estimates

for i=1:N

	y=mu(i,:)';
	b_v=get_biolvar(y,m_groups,biol_var);
	t_v=V(i,:)';
	v=b_v+t_v;

	% null model
	options = gpOptions('ftc');
	options.kern = {'cmpnd','white',{'parametric', struct('variance', {v}, 'input', {x}), 'fixedwhite'},'bias'};
	model = gpCreate(1, 1, x, y, options);
	model1 = gpOptimise(model, 0);
    
	K = model1.K_uu;
	invK = model1.invK_uu;
	logDetK = model1.logDetK_uu;
	ll = -0.5*(logDetK + (model.m)'*invK*(model.m) + size(y, 1)*log(2*pi));
	llLogDet = -.5*(logDetK+size(y, 1)*log(2*pi));
	L0=ll;
     
	% alternative model
	options = gpOptions('ftc');
	options.kern = {'cmpnd','rbf','white',{'parametric', struct('variance', {v}, 'input', {x}), 'fixedwhite'},'bias'};	     
	model = gpCreate(1, 1, x, y, options);
	lengthScale = exp(linspace(log(0.01),log(x(end)),10));
	l=length(lengthScale);
	ll=zeros(l,1);
   	for e = 1:length(lengthScale)
		model.kern.comp{1}.inverseWidth = 1/(lengthScale(e)*lengthScale(e));
		model1 = gpOptimise(model, 0);
		K = model1.K_uu;
		invK = model1.invK_uu;
		logDetK = model1.logDetK_uu;
		ll(e) = -0.5*(logDetK + (model.m)'*invK*(model.m) + size(y, 1)*log(2*pi));
	end
    
	L1=max(ll);
	ff=find(ll==L1);
	ee=ff(1);
	model.kern.comp{1}.inverseWidth = 1/(lengthScale(ee)*lengthScale(ee));
	model1= gpOptimise(model, 0); 
	K = model1.K_uu;
	invK = model1.invK_uu;
	logDetK = model1.logDetK_uu;
	L1 = -0.5*(logDetK + (model.m)'*invK*(model.m) + size(y, 1)*log(2*pi));      
   
    	BF(i)=exp(L1-L0); 
	likli0(i)=exp(L0);
	likli1(i)=exp(L1);
	LS(i)=sqrt(1/(model1.kern.comp{1}.inverseWidth));
    
end

output_fileName='GP_summary';
result=[LS BF likli0 likli1];
dlmwrite(output_fileName, result, 'delimiter','\t')

end
