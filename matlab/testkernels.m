A=1;
L=2;
Sj=1;
Sk=1;
Dj=1;
Dk=1;
sigma2=0.1;
temppoints=[0:5:1280];
kap=zeros(length(temppoints),length(temppoints));
kim=zeros(length(temppoints),length(temppoints));
kam=zeros(length(temppoints),length(temppoints));
for i=1:length(temppoints),
for j=1:length(temppoints),
  kim(i,j)=kernel_instant_m(temppoints(i),temppoints(j),Sj,Sk,Dj,Dk,A,L,sigma2);
  kap(i,j)=kernel_async_p(temppoints(i),temppoints(j),A,L,sigma2);
  kam(i,j)=kernel_async_m(temppoints(i),temppoints(j),Sj,Sk,Dj,Dk,A,L,sigma2);
end;
end;

