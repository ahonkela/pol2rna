function H=compute_hessian(theta,y)
% compute_hessian  calculates the hessian matrix of the posterior likelihood
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

D=-M*psi(1,alpha)+M*psi(1,(alpha+(R/2)));
E=(M/beta)-sum(1./(beta+0.5*C));
F=(-M*alpha/(beta^2))+(alpha+(R/2))*sum(1./((beta+0.5*C).^2));
H=[D E; E F];

end
