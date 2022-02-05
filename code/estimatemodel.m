clear
beep off;
addpath('C:/Users/Valeria/Desktop/CICLO 2020-1/TESIS/Benchmarkmodel/Commodity/functions');
dfolder='C:/Users/Valeria/Desktop/CICLO 2020-1/TESIS/Benchmarkmodel/Commodity/data/';
sfolder='C:/Users/Valeria/Desktop/CICLO 2020-1/TESIS/Benchmarkmodel/Commodity/results/';
clc
maxfile=1; 

%inputs
REPS=5000;%50000;%total reps %6000
BURN=4000;%45000;%burn in
npart=50;
fsize=REPS-BURN;
L=4; %lags of endogenous variables
LH=2; %lags of volatility
LV=2;
CHECK=1;
maxdraws=1000;

T0=120;
irun=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555
file=maxfile;
 sfile=strcat(sfolder,'forecast',num2str(file),'.mat');
 dfile=strcat(dfolder,'dataTERMun12','.xls');

[data0,names]=xlsread(dfile);

data=data0(:,[3 1 2 4]); %interest rate GDP CPI Stock
names=names(:,[3 1 2 4]);
Y=packr(data);
N=cols(data);
EX=(N*(LH))+1;

X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1) ];
Y=Y(L+1:end,:);
X=X(L+1:end,:);









%/////////////////Priors and starting values///////////////////////
y0=Y(1:T0,:);
x0=[];
for j=1:L
x0=[x0 lag0(y0,j)  ];
end
x0=[x0 ones(T0,1)];
y0=y0(L+1:end,:);
x0=x0(L+1:end,:);
beta0=x0\y0; 
e0=y0-x0*beta0;
sigmax0=(e0'*e0)/T0;



%Minnesota type prior for the VAR coefficients via dummy observations
LAMDAP=0.1;
TAUP=0;
EPSILON=1/1000;
EPSILONH=0.001; %closer to 0 tighter prior
RW=0;
[yd,xd,BETA0,SIGMA0]=getdummies(LAMDAP,TAUP,EPSILON,...
    Y,L,EX,LH,EPSILONH,RW);
% BETA0=zeros(N*(N*L+EX),1);
% SIGMA0=eye(N*(N*L+EX))*1000;





%define priors for unconditional mean of transition eq
mubig=zeros(N,1);





c0=chol(sigmax0);
c0=inv(c0./repmat(diag(c0),1,N))'; %A matrix based on training sample
c0(2,1)=-c0(2,1);
% c0(3,1)=-c0(3,1);
p00_a=1;  %prior variances a1



%starting value for the stochastic volatility
y=Y(T0+1:end,:);
x=X(T0+1:end,:);
T=rows(x);

if irun
[outh,outf,outg]=getinitialvol(data,500,400,T0,L);
save initial outh outf outg;
else
    load initial
end
%define priors for variances
vg0=N+1;             %prior df
g0=diag(outg*vg0);  %prior scale paramter for Q[i]~IG(g0,vg0)

h0=outh(2:end,:);


Qbig=diag(outg);



%priors for transition equation
LAMDAPV=0.05;
TAUPV=0;
EPSILONV=100;
EPSILONVY=0.05;
[ydv,xdv,bv0,sv0]=getdummiesVOLF(LAMDAPV,TAUPV,EPSILONV,log(h0),1,LV,EPSILONVY,outf(1:N),outg);


hlast=zeros(T,N*(LH+1));
hlast(:,1:N)=h0;
i=1;
for j=N+1:N:cols(hlast)
hlag=lag0(h0,i);
hlag=packr(hlag);
hlag=[repmat(hlag(1,:),i,1);hlag];
hlast(:,j:j+N-1)=hlag;
i=i+1;
end
hnew=hlast;

xvar=[x(:,1:cols(X)-1) log(hlast(1:end,N+1:end)) x(:,end:end)];
xfix=log(hlast);

hnew=hlast;


varcoef0=xvar\y;
res=y-xvar*varcoef0;
varcoef0=vec(varcoef0)';

mu0=(diag(sigmax0))';
for j=1:LH
mu0=[mu0 (diag(sigmax0))'];  %ln(H0)~N(mu0,s0)
end
s0=eye(cols(mu0))*0.1;
is0=inv(s0);
NS=cols(hlast);

save priors
%Gibbs algorithm
amatsave=zeros(fsize,N,N);
bsave=zeros(fsize,1,N*(N*L+EX));
fsave=zeros(fsize,1,N*(N*LV+(N+1)));
qsave=zeros(fsize,1,N*N);
musave=zeros(fsize,1,N);
hsave=zeros(fsize,T,N);


naccept=zeros(T+1,1);

jgibbs=1;
for igibbs=1:REPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Step 1 A Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%
igibbs
chck=-1;
while chck<0;
amatx=[];
for j=2:N
%v2=-a1*v1+sqrt(exp(h2))*e2
ytemp=res(:,j);  %v2
xtemp=res(:,1:j-1)*(-1); %-v1
ytemp=ytemp./sqrt(hlast(1:rows(hlast),j));
xtemp=xtemp./repmat(sqrt(hlast(1:rows(hlast),j)),1,cols(xtemp));
%prior means and variances
alpha21_0=c0(j,1:j-1)';
% v00_a=diag(abs(alpha21_0))*p00_a;
v00_a=eye(cols(xtemp))*100;
if j<N
v00_a(1,1)=0.001;
end
% alpha21_0
% v00_a
alpha21_2=getreg(ytemp,xtemp,alpha21_0,v00_a,1); %draw from Normal

amatx=[amatx alpha21_2'];
end
amat=repmat(amatx,T,1);
iamat=invpd(chofac(N,amatx'));%chofac converts to a lower triangular 
e1=iamat(2,1)<0; %un
e2=iamat(3,1)<0; %CPI
e3=iamat(4,1)<0; %STOCK

if e1+e2+e3==3;
    chck=1;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Step 2 VAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%coefficients%%%%%%%%%%%%%%%%%%%%%%%%
[beta2,res,roots,problem]=carterkohnvar(y,xvar,0,iamat,[hlast(1,1:N);hlast(:,1:N)],BETA0',SIGMA0,L,CHECK,maxdraws,EX);
if problem
    beta2=varcoef0;
else
    varcoef0=beta2;
    
end

%%%%%%%%%%%%%%%%%%%Step 3 Transition
%%%%%%%%%%%%%%%%%%%equation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  Y0=log(hlast(1:end,1:N));
  
  X0=[lag0(Y0,1) xvar(:,1:N*LV) ones(rows(Y0),1) ];
  Y0=Y0(2:end,:);
  X0=X0(2:end,:);
  Y0s=[Y0;ydv];
  X0s=[X0;xdv];
  iX0=invpd(X0s'*X0s);
  mstar=vec(iX0*(X0s'*Y0s));
  vstar=kron(Qbig,iX0);
  %draw beta but ensure stability
  chck=-1;
  while chck<0
   betaf=mstar+(randn(1,N*(N*LV+(N+1)))*chol(vstar))';
   
   if ~stability(betaf,N,1,N*LV+1)
       chck=10;
   end
  end
  
  
  %step 3 draw MU the long run mean conditional on beta and sigma (see
  %Appendix A in villani.
  betaf1=reshape(betaf,N*LV+(N+1),N);
  
%sample Qbig
resv=Y0-X0*betaf1;
scalev=resv'*resv+g0;
Qbig=iwpq(T+vg0,invpd(scalev));


    

%%%%%%%%%%%%Step 4 Stochastic vol%%%%%%%%
Fmat=zeros(NS,NS);
Fmat(N+1:NS,1:NS-N)=eye(NS-N);
Fmat(1:N,1:N)=betaf1(1:N,:)';
fit=X0(:,N+1:end)*betaf1(N+1:end,:);
fit=[fit(1,:);fit];
Qmat=zeros(NS,NS);
Qmat(1:N,1:N)=Qbig;
iQmat=zeros(NS,NS);
iQmat(1:N,1:N)=invpd(Qbig);
cQmat=zeros(NS,NS);
cQmat(1:N,1:N)=chol(Qbig);
varcoef=reshape(beta2,N*L+EX,N); 

%Sample volatility
[betax,BB,W]=pftest3F(Fmat,fit,cQmat,iQmat,iamat,xfix,y,x,npart,NS,log(mu0),varcoef,N);
index=discretesample(W(:,end),1);
xfix=betax(:,:,index);

%update
hlast=exp(xfix);
xvar=[x(:,1:cols(X)-1) log(hlast(1:end,N+1:end)) x(:,end:end)];

%impulse response to interest rate volatility
if igibbs>BURN
  
    
      amatsave(jgibbs,:,:)=iamat;
bsave(jgibbs,:,:)=beta2;
fsave(jgibbs,:,:)=betaf';
qsave(jgibbs,:,:)=vec(Qbig)';
musave(jgibbs,:,:)=mubig';
hsave(jgibbs,:,:)=hlast(:,1:N);
    
    
    
    jgibbs=jgibbs+1;
end

disp(sprintf('Iteration Number= %s file no=%s ', num2str([igibbs min(naccept) median(naccept)]),num2str(file)));
if rem(igibbs,100)==0
    save progress;
end
end
 save(sfile,'amatsave','bsave','fsave','qsave','musave','hsave','naccept'); 
 
getirfmcnew