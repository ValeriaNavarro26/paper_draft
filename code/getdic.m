clear
beep off;
addpath('G:\Mi unidad\CICLO 2020-1\TESIS\Benchmarkmodel\Commodity/functions');
dfolder='G:\Mi unidad\CICLO 2020-1\TESIS\Benchmarkmodel\Commodity/data/';
sfolder='G:\Mi unidad\CICLO 2020-1\TESIS\Benchmarkmodel\Commodity//results/';
clc
file=1; 

horizon=60;
load priors
 sfile=strcat(sfolder,'forecast',num2str(file),'.mat');
load(sfile);
npartx=500; %number of particles
fsize=500; %number of reps




lik=zeros(fsize,1);

parfor jgibbs=1:fsize
    jgibbs
    iamat=squeeze(amatsave(jgibbs,:,:));
    beta2=bsave(jgibbs,:,:);
    varcoef=reshape(beta2,N*L+EX,N);
fbig=squeeze(fsave(jgibbs,:,:));
Qbig=squeeze(qsave(jgibbs,:,:));
mubig=squeeze(musave(jgibbs,:,:));
hlast=squeeze(hsave(jgibbs,:,:));
betaf1=reshape(fbig,N*LV+(N+1),N);

B00=zeros(1,NS);
kk=1;
for k=1:LH+1
    B00(1,kk:kk+N-1)=log(hlast(LH-k+2,:));
    kk=kk+N;
end
P00=eye(cols(B00)).*0.1;

Fmat=zeros(NS,NS);
Fmat(N+1:NS,1:NS-N)=eye(NS-N);
Fmat(1:N,1:N)=betaf1(1:N,:)';
X0=[xvar(:,1:N*LV) ones(rows(xvar),1)];
fit=X0*betaf1(N+1:end,:);
Qmat=zeros(NS,NS);
Qmat(1:N,1:N)=chol(reshape(Qbig,N,N));






[ likx,states ] = ...
    particlefilterlinear( y,x,varcoef,iamat,Fmat,...
    Qmat,fit,...
    npartx,rows(y),N,L,EX,B00,P00);
lik(jgibbs)=likx;

jgibbs


 
end

%posterior mean
iamatm=squeeze(mean(amatsave(1:fsize,:,:)));
    beta2m=squeeze(mean(bsave(1:fsize,:,:)));
    varcoefm=reshape(beta2m,N*L+EX,N);
fbigm=squeeze(mean(fsave(1:fsize,:,:)));
Qbigm=squeeze(mean(qsave(1:fsize,:,:)));
mubigm=squeeze(mean(musave(1:fsize,:,:)));
hlastm=squeeze(mean(hsave(1:fsize,:,:)));
betaf1=reshape(fbigm,N*LV+(N+1),N);

B00=zeros(1,NS);
kk=1;
for k=1:LH+1
    B00(1,kk:kk+N-1)=log(hlastm(LH-k+2,:));
    kk=kk+N;
end
P00=eye(cols(B00)).*0.1;

Fmat=zeros(NS,NS);
Fmat(N+1:NS,1:NS-N)=eye(NS-N);
Fmat(1:N,1:N)=betaf1(1:N,:)';
X0=[xvar(:,1:N*LV) ones(rows(xvar),1)];
fit=X0*betaf1(N+1:end,:);
Qmat=zeros(NS,NS);
Qmat(1:N,1:N)=chol(reshape(Qbigm,N,N));










[ likm,statesm ] = ...
    particlefilterlinear( y,x,varcoefm,iamatm,Fmat,...
    Qmat,fit,...
    npartx,rows(y),N,L,EX,B00,P00);

params=mean(-2*lik)-(-2*likm);
   bdic1=mean(-2*lik)+params;

save('bdic1','bdic1','params')
clear all

