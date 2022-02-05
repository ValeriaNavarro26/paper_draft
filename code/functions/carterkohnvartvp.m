function[beta2,error,roots,epsilon,problem]=...
    carterkohnvartvp(Y,X,Q,iamat,hlast,lamda,beta0,P00,L,CHECK,maxdraws,EX)
T=rows(Y);
N=cols(Y);
amat=invpd(iamat);
%%Step 2a Set up matrices for the Kalman Filter

ns=cols(beta0);
F=eye(ns);
mu=0;
beta_tt=zeros(T,ns);          %will hold the filtered state variable
ptt=zeros(T,ns,ns);    % will hold its variance
beta11=beta0;
p11=P00;


% %%%%%%%%%%%Step 2b run Kalman Filter

for i=1:T
   x=kron(eye(N),X(i,:));

H=diag(hlast(i+1,:).*(1./lamda(i,:)));
R=iamat*H*iamat';
    %Prediction
beta10=beta11;
p10=p11+Q;
yhat=(x*(beta10)')';                                               
eta=Y(i,:)-yhat;
feta=(x*p10*x')+R;
%updating
K=(p10*x')*invpd(feta);
beta11=(beta10'+K*eta')';
p11=p10-K*(x*p10);

ptt(i,:,:)=p11;
beta_tt(i,:)=beta11;

end


%%%%%%%%%%%end of Kalman Filter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%step 2c Backward recursion to calculate the mean and variance of the distribution of the state
%vector
chck=-1;
problem=0;
trys=1;
while chck<0 && trys<=maxdraws
    
beta2 = zeros(T,ns);   %this will hold the draw of the state variable
wa=randn(T,ns);
error=zeros(T,N);
epsilon=zeros(T,N);
roots=zeros(T,1);

i=T;  %period t
p00=squeeze(ptt(i,:,:)); 
beta2(i,:)=beta_tt(i:i,:)+(wa(i:i,:)*cholx(p00));   %draw for beta in period t from N(beta_tt,ptt)
error(i,:)=Y(i,:)-X(i,:)*reshape(beta2(i:i,:),N*L+EX,N);  %var residuals
epsilon(i,:)=error(i,:)*amat';
roots(i)=stability(beta2(i,:)',N,L,EX);

%periods t-1..to .1

for i=T-1:-1:1
pt=squeeze(ptt(i,:,:));
iFptF=invpd(pt+Q);
% iFptFx=invpd(F*pt*F'+Q);
% iFptF-iFptFx
bm=beta_tt(i:i,:)+(pt*iFptF*(beta2(i+1:i+1,:)-beta_tt(i,:))')';  %update the filtered beta for information contained in beta[t+1]    
% bmx=beta_tt(i:i,:)+(pt*F'*iFptF*(beta2(i+1:i+1,:)-beta_tt(i,:)*F')')';  %update the filtered beta for information contained in beta[t+1]                                                                                 %i.e. beta2(i+1:i+1,:) eq 8.16 pp193 in Kim Nelson
% bm-bmx

%i.e. beta2(i+1:i+1,:) eq 8.16 pp193 in Kim Nelson
pm=pt-pt*iFptF*pt;  %update covariance of beta
% pmx=pt-pt*F'*iFptF*F*pt;  %update covariance of beta
% pm-pmx

beta2(i:i,:)=bm+(wa(i:i,:)*cholx(pm));  %draw for beta in period t from N(bm,pm)eq 8.17 pp193 in Kim Nelson
error(i,:)=Y(i,:)-X(i,:)*reshape(beta2(i:i,:),N*L+EX,N);  %var residuals
epsilon(i,:)=error(i,:)*amat';
roots(i)=stability(beta2(i,:)',N,L,EX);
end
if CHECK
if sum(roots)==0
    chck=1;
else
    trys=trys+1;
end
else
 chck=1;
end
end
if CHECK
    if chck<0
        problem=1;
    end
end