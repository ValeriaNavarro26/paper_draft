function [yhatx,hhatx,yhatz,hhatz,yhat,hhat]=getdatmp(TT,hlast,y,L,LV,Qbig,betaf1,varcoef,iamat,LH)
N=cols(y);
yhatx=zeros(TT,N);
yhatx(1:L,:)=y(1:L,:);
hhatx=zeros(TT,N);
hhatx(1:L,:)=log(hlast(1:L,:));
% all level shocks
yhatz=yhatx;
hhatz=hhatx;
%benchmark
yhat=yhatx;
hhat=hhatx;


for t=L+1:TT
    %SVOL
    residH=randn(1,N)*chol(reshape(Qbig,N,N));
    xh=hhatx(t-1,:);
    xz=hhatz(t-1,:);
    x=hhat(t-1,:);
    for j=1:LV
        xh=[xh yhatx(t-j,:)];
        xz=[xz yhatz(t-j,:)];
        x=[x yhat(t-j,:)];
    end
            xh=[xh 1];
            xz=[xz 1];
            x=[x 1];
            hhatx(t,:)=xh*betaf1;
            hh=exp(hhatx(t,:));
            hhatz(t,:)=xz*betaf1;
            hhz=exp(hhatz(t,:));
            hhat(t,:)=x*betaf1+residH;
            h=exp(hhat(t,:));
    %structural shocks
    A0=iamat*diag(sqrt(hh));
    sres=randn(1,N);
    A0z=iamat*diag(sqrt(hhz));
    a0=iamat*diag(sqrt(h));
    %shut down all shocks except MP
    sresx=sres;
    sresx(2:4)=0;
    residx=sresx*A0';
    
    residz=sres*A0z';
    
    resid=sres*a0';
   
    xhatx=[];
    xhatz=[];
    xhat=[];
    for j=1:L
        xhatx=[xhatx yhatx(t-j,:)];
        xhatz=[xhatz yhatz(t-j,:)];
        xhat=[xhat yhat(t-j,:)];
    end
    for j=1:LH
        xhatx=[xhatx hhatx(t-j,:)];
        xhatz=[xhatz hhatz(t-j,:)];
        xhat=[xhat hhat(t-j,:)];
    end
    xhatx=[xhatx 1];
    xhatz=[xhatz 1];
    xhat=[xhat 1];
    
    yhatx(t,:)=xhatx*varcoef+residx;
    yhatz(t,:)=xhatz*varcoef+residz;
    yhat(t,:)=xhat*varcoef+resid;
    
    
end