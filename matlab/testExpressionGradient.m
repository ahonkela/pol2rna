delay=disimKern.delay;
iWidth=disimKern.inverseWidth;
t1=jointmodel.t;
t2=jointmodel.t;
dim1 = size(t1, 1);
dim2 = size(t2, 1);
sigma = sqrt(2/iWidth);
sigmaA=sqrt(2/disimKern.inverseWidth);
D_i = disimKern.decay;
variancemultiplier=disimKern.di_variance*sqrt(disimKern.variance);
origt1=t1;origt1PosFlag=origt1>delay;origt1PosFlag=origt1PosFlag(:, ones(1, dim2));
t1=t1-delay;
I=find(t1<0);
t1(I)=0;
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

%value1=-sqrt(pi)*sigma/(2*D_i)*...
%    ((diffT-1/D_i).*erf(diffT/sigma)
%    -(t1Mat-1/D_i).*erf(t1Mat/sigma)
%    -(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma));

tempgrad=...
    +sqrt(pi)*sigma/(2*(D_i^2))*...
    ( (diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma) ) ...
    -sqrt(pi)*sigma/(2*D_i)*...
    (-1/(D_i^2)*exp(lnDiffErfs(t1Mat/sigma,diffT/sigma)) ...
     +((t1Mat+1/D_i).*exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma) );
    
%    +sqrt(pi)*sigma/(2*(D_i^2))*...
%    ((diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma));
%    -sqrt(pi)*sigma/(2*D_i)*...
%    (-1/(D_i^2)*exp(lnDiffErfs(t1Mat/sigma,diffT/sigma)) 
%     +((t1Mat+1/D_i).*exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma) );


%value2=-sqrt(pi)*sigma/(2*D_i)*...
%    ((diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma));

[lndiff1,lndiffsigns1]=lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2);
dK_dD_part1a=...
   +2*sqrt(pi)*sigma/(2*(D_i^3))*...
    ( lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lndiff1)...
     +lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) ) ...
   -sqrt(pi)*sigma/(2*D_i*D_i)*...
    ( lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lndiff1) .*(D_i*(sigma^2)/2 -t1Mat+t2Mat) ...
     +( exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat -(D_i*sigma/2+t2Mat/sigma).^2) ...
       -exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat -(D_i*sigma/2+t2Mat/sigma-t1Mat/sigma).^2) )*sigma/sqrt(pi) ...
     ...
     +lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) .*(D_i*(sigma^2)/2 -t1Mat) ...
     + (  exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2-t1Mat/sigma).^2) ...
        - exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2).^2) )*sigma/sqrt(pi)   );


tempgrad=...
   +2*sqrt(pi)*sigma/(2*(D_i^3))*...
    ( exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma))...
     +exp((D_i*sigma/2)^2-D_i*t1Mat+lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2)) ) ...
   -sqrt(pi)*sigma/(2*D_i*D_i)*...
    ( exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma))...
       .*(D_i*(sigma^2)/2 -t1Mat+t2Mat) ...
     +( exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat -(D_i*sigma/2+t2Mat/sigma).^2) ...
       -exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat -(D_i*sigma/2+t2Mat/sigma-t1Mat/sigma).^2) )*sigma/sqrt(pi) ...
     ...
     +exp((D_i*sigma/2)^2-D_i*t1Mat+lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2)) .*(D_i*(sigma^2)/2 -t1Mat) ...
     + (  exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2-t1Mat/sigma).^2) ...
        - exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2).^2) )*sigma/sqrt(pi)   );
    

tempgrad=...
    -sqrt(pi)*sigma/(2*D_i)*origt1PosFlag.* ...
    ( exp(lnDiffErfs(t1Mat/sigma,diffT/sigma)) -exp(-D_i*t1Mat).*erf(t2Mat/sigma) ) ...
    -1/D_i*origt1PosFlag.* ...
    ( -(diffT-1/D_i).*exp(-(diffT/sigma).^2) +(t1Mat-1/D_i).*exp(-(t1Mat/sigma).^2) );

    -sqrt(pi)*sigma/(2*D_i)*origt1PosFlag.* ...
    (-erf(diffT/sigma) -(diffT-1/D_i).*exp(-(diffT/sigma).^2)/sigma*2/sqrt(pi) ...
     +erf(t1Mat/sigma) +(t1Mat-1/D_i).*exp(-(t1Mat/sigma).^2)/sigma*2/sqrt(pi) ...
     -exp(-D_i*t1Mat).*erf(t2Mat/sigma)  );
    
    -sqrt(pi)*sigma/(2*D_i)*origt1PosFlag.* ...
    ( exp(lnDiffErfs(t1Mat/sigma,diffT/sigma)) -exp(-D_i*t1Mat).*erf(t2Mat/sigma) ) ...
    -1/D_i*origt1PosFlag.* ...
    ( -(diffT-1/D_i).*exp(-(diffT/sigma).^2) +(t1Mat-1/D_i).*exp(-(t1Mat/sigma).^2) );


    
     -1/D_i*origt1PosFlag .* ( exp(-(diffT/sigma).^2).*diffT - exp(-(t1Mat/sigma).^2).*t1Mat );
    
    
    -sqrt(pi)*sigma/(2*D_i*D_i)*origt1PosFlag.*...
    ( exp(((D_i*sigma/2)^2)-D_i*(t1Mat-t2Mat)   +lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma))*D_i ...
      -2/(sigma*sqrt(pi))*exp((D_i*sigma/2)^2-D_i*(t1Mat-t2Mat) -(D_i*sigma/2+(t2Mat-t1Mat)/sigma).^2) ...
     ...
     +exp((D_i*sigma/2)^2-D_i*t1Mat+lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2))*D_i ...
     +2/(sigma*sqrt(pi))*exp((D_i*sigma/2)^2-D_i*t1Mat -(D_i*sigma/2-t1Mat/sigma).^2)  );
    

[lndiff1,lndiffsigns1]=lnDiffErfs(t1Mat/sigma,diffT/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(t2Mat/sigma,-diffT/sigma);
tempgrad=...
    -sqrt(pi)/(2*D_i)*...
    (-t1Mat.*lndiffsigns1.*exp(lndiff1) -t2Mat.*lndiffsigns2.*exp(lndiff2) ...
     +1/D_i*(erf(t1Mat/sigma)-erf(diffT/sigma)) - (exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma) ) ...
    -1/(D_i*sigma)*...
    ( (diffT-1/D_i).*(-diffT).*exp(-(diffT/sigma).^2) ...
     -(t1Mat-1/D_i).*(-t1Mat).*exp(-(t1Mat/sigma).^2) ...
     -(t2Mat+exp(-D_i*t1Mat)/D_i).*(-t2Mat).*exp(-(t2Mat/sigma).^2) );

tempgrad=tempgrad*(-disimKern.inverseWidth^(-1.5)/sqrt(2));


tempgrad=...
    -sigma/D_i*...
    (exp(-(diffT/sigma).^2)-exp(-(t1Mat/sigma).^2)+1-exp(-(t2Mat/sigma).^2)) ...
    -1/(sigma*D_i)*...
    ((diffT.^2).*exp(-(diffT/sigma).^2)-(t1Mat.^2).*exp(-(t1Mat/sigma).^2)-(t2Mat.^2).*exp(-(t2Mat/sigma).^2));


% -sqrt(pi)*sigma/(2*D_i*D_i)*...
%    (exp((D_i*sigma/2)^2-D_i*t1Mat+D_i*t2Mat+lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+t2Mat/sigma-t1Mat/sigma)) 
%     +exp((D_i*sigma/2)^2-D_i*t1Mat+lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2)) );


[lndiff1,lndiffsigns1]=lnDiffErfs(D_i*sigma/2+t2Mat/sigma,D_i*sigma/2+(t2Mat-t1Mat)/sigma);
[lndiff2,lndiffsigns2]=lnDiffErfs(D_i*sigma/2-t1Mat/sigma,D_i*sigma/2);
tempgrad=...
    -(sqrt(pi)/2)*(1/(D_i*D_i) + (sigma^2)/2)*...
     lndiffsigns1.*exp((D_i*sigma/2)^2-D_i*(t1Mat-t2Mat)+lndiff1) ...
    -1/D_i * exp((D_i*sigma/2)^2 - D_i*(t1Mat-t2Mat)) .* ...
      (  (sigma/2-t2Mat/(sigma*D_i)).*exp(-(D_i*sigma/2+t2Mat/sigma).^2) ...
       - (sigma/2+(t1Mat-t2Mat)/(sigma*D_i)).*exp(-(D_i*sigma/2+(t2Mat-t1Mat)/sigma).^2) ) ...
    ...
    -(sqrt(pi)/2)*(1/(D_i*D_i) + (sigma^2)/2)*...
     lndiffsigns2.*exp((D_i*sigma/2)^2-D_i*t1Mat+lndiff2) ...
    -1/D_i * exp((D_i*sigma/2)^2 - D_i*t1Mat).* ...
     (   (sigma/2+t1Mat/(sigma*D_i)).*exp(-(D_i*sigma/2-t1Mat/sigma).^2) ...
       - (sigma/2)*exp(-(D_i*sigma/2)^2) ) ;

tempgrad=tempgrad*(-disimKern.inverseWidth^(-1.5)/sqrt(2));


myeps=1e-5;

D_i=disimKern.decay-myeps;
D_iA=disimKern.decay;
iWidth=disimKern.inverseWidth;
delay=disinKern.delay;
t1=jointmodel.t;
t2=jointmodel.t;
dim1 = size(t1, 1);
dim2 = size(t2, 1);
sigma = sqrt(2/iWidth);
sigmaA=sqrt(2/disimKern.inverseWidth);
variancemultiplier=disimKern.di_variance*sqrt(disimKern.variance);
origt1=t1;origt1PosFlag=origt1>delay;origt1PosFlag=origt1PosFlag(:, ones(1, dim2));
t1=t1-delay;
I=find(t1<0);
t1(I)=0;
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
value1=-sqrt(pi)*sigma/(2*D_i)*...
    ((diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma));


D_i=disimKern.decay+myeps;
D_iA=disimKern.decay;
iWidth=disimKern.inverseWidth;
t1=jointmodel.t;
t2=jointmodel.t;
dim1 = size(t1, 1);
dim2 = size(t2, 1);
sigma = sqrt(2/iWidth);
sigmaA=sqrt(2/disimKern.inverseWidth);
delay = disimKern.delay;
variancemultiplier=disimKern.di_variance*sqrt(disimKern.variance);
origt1=t1;origt1PosFlag=origt1>delay;origt1PosFlag=origt1PosFlag(:, ones(1, dim2));
t1=t1-delay;
I=find(t1<0);
t1(I)=0;
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
value2=-sqrt(pi)*sigma/(2*D_i)*...
    ((diffT-1/D_i).*erf(diffT/sigma)-(t1Mat-1/D_i).*erf(t1Mat/sigma)-(t2Mat+exp(-D_i*t1Mat)/D_i).*erf(t2Mat/sigma));

diffgrad=real(value2-value1)/(2*myeps);




























t1=jointmodel.t;
t2=jointmodel.t;

sigma = sqrt(2/kern.inverseWidth);
if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    variancemultiplier = (kern.sensitivity*kern.sensitivity);
else
    variancemultiplier = kern.variance;
end

dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

k1a = sqrt(pi)/2*( t1Mat.*erf(t1Mat/sigma) + t2Mat.*erf(t2Mat/sigma) - diffT.*erf(diffT/sigma) );
k1b = sigma/2*( exp(-(t1Mat/sigma).^2) + exp(-(t2Mat/sigma).^2) - exp(-(diffT/sigma).^2) - 1 );
K = (k1a+k1b)*sigma*variancemultiplier;

deriv_variancemultiplier = sum(sum(sigma*(k1a+k1b).*covGrad));

%deriv_sigma = (k1a+k1b)*variancemultiplier ...
%    + variancemultiplier*sigma/2*(exp(-(t1Mat/sigma).^2)+exp(-(t2Mat/sigma).^2)-exp(-(diffT/sigma).^2)-1) ...
%    + variancemultiplier/sigma*(-(t1Mat.^2).*exp(-(t1Mat/sigma).^2)-(t2Mat.^2).*exp(-(t2Mat/sigma).^2)+(diffT.^2).*exp(-(diffT/sigma).^2)) ...
%    + variancemultiplier/sigma*(-(t1Mat.^2).*exp(-(t1Mat/sigma).^2)-(t2Mat.^2).*exp(-(t2Mat/sigma).^2)+(diffT.^2).*exp(-(diffT/sigma).^2)) ;

deriv_sigma = (k1a+2*k1b)*variancemultiplier;% ...
%    + variancemultiplier/sigma*(-(t1Mat.^2).*exp(-(t1Mat/sigma).^2)-(t2Mat.^2).*exp(-(t2Mat/sigma).^2)+(diffT.^2).*exp(-(diffT/sigma).^2)) ...
%    + variancemultiplier/sigma*(+(t1Mat.^2).*exp(-(t1Mat/sigma).^2)+(t2Mat.^2).*exp(-(t2Mat/sigma).^2)-(diffT.^2).*exp(-(diffT/sigma).^2)) ;

%deriv_sigma = (2*k1b)*variancemultiplier ...
%    + variancemultiplier/(sigma)*(+(t1Mat.^2).*exp(-(t1Mat/sigma).^2)+(t2Mat.^2).*exp(-(t2Mat/sigma).^2)-(diffT.^2).*exp(-(diffT/sigma).^2)) ;



%    + variancemultiplier*sigma/2*(exp(-(t1Mat/sigma).^2)+exp(-(t2Mat/sigma).^2)-exp(-(diffT/sigma).^2)-1) ...

deriv_inversewidth=deriv_sigma*(-0.5*sigma/kern.inverseWidth);
%deriv_sigma = sum(sum(deriv_sigma.*covGrad));





myeps=1e-3;


iWidth=kern.inverseWidth-myeps;
sigma = sqrt(2/iWidth);
if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    variancemultiplier = (kern.sensitivity*kern.sensitivity);
else
    variancemultiplier = kern.variance;
end
t1=jointmodel.t;
t2=jointmodel.t;
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
k1a = sqrt(pi)/2*( t1Mat.*erf(t1Mat/sigma) + t2Mat.*erf(t2Mat/sigma) - diffT.*erf(diffT/sigma) );
k1b = sigma/2*( exp(-(t1Mat/sigma).^2) + exp(-(t2Mat/sigma).^2) - exp(-(diffT/sigma).^2) - 1 );
value1 = (k1a+k1b)*sigma*variancemultiplier;
%value1 = (k1b)*sigma*variancemultiplier;

iWidth=kern.inverseWidth+myeps;
sigma = sqrt(2/iWidth);
if isfield(kern, 'isNegativeS') && (kern.isNegativeS == true)
    variancemultiplier = (kern.sensitivity*kern.sensitivity);
else
    variancemultiplier = kern.variance;
end
t1=jointmodel.t;
t2=jointmodel.t;
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
k1a = sqrt(pi)/2*( t1Mat.*erf(t1Mat/sigma) + t2Mat.*erf(t2Mat/sigma) - diffT.*erf(diffT/sigma) );
k1b = sigma/2*( exp(-(t1Mat/sigma).^2) + exp(-(t2Mat/sigma).^2) - exp(-(diffT/sigma).^2) - 1 );
value2 = (k1a+k1b)*sigma*variancemultiplier;
%value2 = (k1b)*sigma*variancemultiplier;

diffgrad=(value2-value1)/(2*myeps)







D=kern.decay;
l=sqrt(2/kern.inverseWidth);
delay=kern.delay;
origt1=jointmodel.t;t1=origt1-delay;I=find(t1<0);t1(I)=0;
origt2=jointmodel.t;t2=origt2-delay;I=find(t2<0);t2(I)=0;

dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);
origt1PosFlag=origt1>delay;origt1PosFlag=origt1PosFlag(:, ones(1, dim2));
origt2PosFlag=origt2>delay;origt2PosFlag=origt2PosFlag(:, ones(1, dim1))';

%tempgrad=...
%  0.5*sqrt(pi)/(D^2) * ...
%    (  -l*origt1PosFlag.*erf(t1Mat/l)  + l*origt2PosFlag.*exp(-D*t2Mat +log(erf(t1Mat/l)))  ...
%       -origt1PosFlag.*(t1Mat - 1/D).*exp(-(t1Mat/l).^2) 
%       -origt1PosFlag.*exp(-log(D) -D*t2Mat -(t1Mat/l).^2) );
    
tempgrad=...   
    sqrt(pi)*l * ...
     0.25*exp(-2*log(D) + (D*l/2)^2 - D*(diffT) + lnDiffErfs(D*l/2+t2Mat/l,D*l/2-diffT/l)).*(-origt2PosFlag+origt1PosFlag) ...
  + 0.5*exp(-3*log(D) + (D*l/2)^2 - D*(diffT) ) ...
    .*( exp(-(D*l/2+t2Mat/l).^2).*(-origt2PosFlag) -exp(-(D*l/2-diffT/l).^2).*(-origt2PosFlag+origt1PosFlag) ) ...
    + sqrt(pi)*l * ...
     0.25*exp(-2*log(D) + (D*l/2)^2 - D*(t1Mat+t2Mat) + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)).*(origt2PosFlag+origt1PosFlag) ...
  + 0.5*exp(-3*log(D) + (D*l/2)^2 - D*(t1Mat+t2Mat) ) ...
    .*(  -exp(-(D*l/2-t1Mat/l).^2).*(origt1PosFlag) ) ...
    - sqrt(pi)*l * ...
     0.5*exp(-2*log(D) + (D*l/2)^2 - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)).*(origt1PosFlag) ...
  - exp(-3*log(D) + (D*l/2)^2 - D*(t1Mat) ) ...
    .*(  -exp(-(D*l/2-t1Mat/l).^2).*(origt1PosFlag) ) ;



    sqrt(pi)*l * ...
    ( 0.25*exp(-3*log(D) + (D*l/2)^2 + D*(diffT) + lnDiffErfs(D*l/2+t1Mat/l,D*l/2+diffT/l)) ...
     +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
     -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) );







    
    
  0.25*sqrt(pi)*l*kern.di_variance*kern.variance * ...
    ( -3*exp(-4*log(D) + (D*l/2)^2 + D*(diffT) + lnDiffErfs(D*l/2+t1Mat/l,D*l/2+diffT/l)) ...
      + 0.5*(l^2)*exp(-2*log(D) + (D*l/2)^2 + D*(diffT) + lnDiffErfs(D*l/2+t1Mat/l,D*l/2+diffT/l)) ...
      + diffT.*exp(-3*log(D) + (D*l/2)^2 + D*(diffT) + lnDiffErfs(D*l/2+t1Mat/l,D*l/2+diffT/l)) ...
      + exp(-3*log(D) + (D*l/2)^2 + D*(diffT) -(D*l/2+t1Mat/l).^2)*(l/sqrt(pi)) ...
      - exp(-3*log(D) + (D*l/2)^2 + D*(diffT) -(D*l/2+diffT/l).^2)*(l/sqrt(pi)) ...
      ...
      -3*exp(-4*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      +0.5*(l^2)*exp(-2*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      -(t2Mat+t1Mat).*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      +exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat -(D*l/2)^2)*(l/sqrt(pi)) ...
      -exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat -(D*l/2-t2Mat/l).^2)*(l/sqrt(pi)) ...
      ...      
      +6*exp(-4*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...           
      -(l^2)*exp(-2*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      +2*t2Mat.*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      -2*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - (D*l/2)^2)*(l/sqrt(pi)) ...
      +2*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - (D*l/2-t2Mat/l).^2)*(l/sqrt(pi)) );



      ( 0.25*exp(-3*log(D) + (D*l/2)^2 + D*(diffT) + lnDiffErfs(D*l/2+t1Mat/l,D*l/2+diffT/l)) ...
       +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
      -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) );





    
    


%      +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ...
%      -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t2Mat/l)) ) ...
%      +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat - D*t1Mat).* ...
%	((2/sqrt(pi))*exp(-(D*l/2)^2)*(D/2)-(2/sqrt(pi))*exp(-(D*l/2-t2Mat/l).^2).*(D/2+t2Mat/(l^2))) ...
%      -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t2Mat).* ...
%	((2/sqrt(pi))*exp(-(D*l/2)^2)*(D/2)-(2/sqrt(pi))*exp(-(D*l/2-t2Mat/l).^2).*(D/2+t2Mat/(l^2))) );
    
    

myeps=1e-1;
D=kern.decay;
l=sqrt(2/kern.inverseWidth);
delay0=kern.delay;

delay=delay0-myeps;
t2=jointmodel.t-delay;I=find(t2<0);t2(I)=0;
t1=jointmodel.t-delay;I=find(t1<0);t1(I)=0;
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

value1=... 
sqrt(pi)*l * ...
    ( 0.25*exp(-3*log(D) + (D*l/2)^2 - D*(diffT) + lnDiffErfs(D*l/2+t2Mat/l,D*l/2-diffT/l)) ...
     +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t1Mat - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)) ...
     -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)) );    

delay=delay0+myeps;
t2=jointmodel.t-delay;I=find(t2<0);t2(I)=0;
t1=jointmodel.t-delay;I=find(t1<0);t1(I)=0;
dim1 = size(t1, 1);
dim2 = size(t2, 1);
t1Mat = t1(:, ones(1, dim2));
t2Mat = t2(:, ones(1, dim1))';
diffT = (t1Mat - t2Mat);

value2=...      
sqrt(pi)*l * ...
    ( 0.25*exp(-3*log(D) + (D*l/2)^2 - D*(diffT) + lnDiffErfs(D*l/2+t2Mat/l,D*l/2-diffT/l)) ...
     +0.25*exp(-3*log(D) + (D*l/2)^2 - D*t1Mat - D*t2Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)) ...
     -0.5*exp(-3*log(D) + (D*l/2)^2 - D*t1Mat + lnDiffErfs(D*l/2,D*l/2-t1Mat/l)) );    

diffgrad=(value2-value1)/(2*myeps)





