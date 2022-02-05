function [out,beta_tt]=getfactor(data,z,H,R,F1,F2,beta0,P00,hlast,iamat1,iamat2,Sbig,exo1,exo2,N,nfact,e1)
T=rows(data);
ns=cols(P00);


%Carter and Kohn algorithm to draw the factor
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
% %%%%%%%%%%%Step 6a run Kalman Filter
i=1;
x=H;
if e1(i)==1
QQ=iamat1*diag(hlast(i+1,1).*Sbig)*iamat1';
Q=zeros(ns,ns);
Q(1:N,1:N)=QQ;
MUx=zeros(1,ns);
MUx(:,1:N)=exo1(i,:);
F=F1;
else
QQ=iamat2*diag(hlast(i+1,1).*Sbig)*iamat2';
Q=zeros(ns,ns);
Q(1:N,1:N)=QQ;
MUx=zeros(1,ns);
MUx(:,1:N)=exo2(i,:);
F=F2;  
end
%Prediction
beta10=MUx+beta0*F';
p10=F*P00*F'+Q;
yhat=(x*(beta10)')';                                                
eta=[data(i,:) z(i,1)]-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
beta_tt(i,:)=beta11;
ptt(i,:,:)=p11;
for i=2:T
    %
  
if e1(i)==1
QQ=iamat1*diag(hlast(i+1,1).*Sbig)*iamat1';
Q=zeros(ns,ns);
Q(1:N,1:N)=QQ;
MUx=zeros(1,ns);
MUx(:,1:N)=exo1(i,:);
F=F1;
else
QQ=iamat2*diag(hlast(i+1,1).*Sbig)*iamat2';
Q=zeros(ns,ns);
Q(1:N,1:N)=QQ;
MUx=zeros(1,ns);
MUx(:,1:N)=exo2(i,:);
F=F2;  
end
    %Prediction
beta10=MUx+beta11*F';
p10=F*p11*F'+Q;
yhat=(x*(beta10)')';                                               
eta=[data(i,:) z(i,1)]-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);
ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;
end
% Backward recursion to calculate the mean and variance of the distribution of the state
%vector
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
jv=1:nfact; %index of state variables to extract
wa=randn(T,ns);

i=T;  %period t
p00=squeeze(ptt(i,jv,jv)); 
beta2(i,jv)=beta_tt(i:i,jv)+(wa(i:i,jv)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
%periods t-1..to .1
for i=T-1:-1:1
    
if e1(i+1)==1
QQ=iamat1*diag(hlast(i+2,1).*Sbig)*iamat1';
Q=zeros(ns,ns);
Q(1:N,1:N)=QQ;
MUx=zeros(1,ns);
MUx(:,1:N)=exo1(i+1,:);
F=F1;
else
QQ=iamat2*diag(hlast(i+2,1).*Sbig)*iamat2';
Q=zeros(ns,ns);
Q(1:N,1:N)=QQ;
MUx=zeros(1,ns);
MUx(:,1:N)=exo2(i+1,:);
F=F2;  
end

q=Q(jv,jv);
mu=MUx(jv);
f=F(jv,:);
pt=squeeze(ptt(i,:,:));
bm=beta_tt(i:i,:)+(pt*f'*invpd(f*pt*f'+q)*(beta2(i+1:i+1,jv)-mu-beta_tt(i,:)*f')')';  
pm=pt-pt*f'*invpd(f*pt*f'+q)*f*pt;  
beta2(i:i,jv)=bm(jv)+(wa(i:i,jv)*cholx(pm(jv,jv)));  

end

out=beta2(:,jv);