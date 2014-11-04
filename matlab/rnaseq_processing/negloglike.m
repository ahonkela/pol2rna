function negloglik=negloglike(theta,y)
% theta = [alpha; beta]
% y = the data vector, M*R matrix for M transcripts R replicates

n = size(y);
M=n(1);
R=n(2);
alpha=theta(1);
beta=theta(2);
mu=mean(y')';
y_2=y.^2;
A=sum(y_2')';
B=R*mu.^2;
C=A-B;

negloglik=-(M*alpha*log(beta)-M*gammaln(alpha)+M*gammaln(alpha+R/2)-(alpha+R/2)*sum(log(beta+0.5*C)));

end
