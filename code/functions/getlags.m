function L=getlags(data,maxlag,choice)

N=cols(data);

%select lag length
sic=zeros(maxlag,1);
aic=zeros(maxlag,1);
for i=1:maxlag
    Y=data;
    L=i;
    
    p=N*L+1;
    NN=N*(N*L+1);
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(rows(X),1)];
Y=Y(L+1:end,:);
X=X(L+1:end,:);
T=rows(X);
beta2=X\Y;
res=Y-X*beta2;
sigma=(res'*res)/(T-p);
detsigma=det(sigma);
lik=-T/2*(N*(1+log(2*pi))+log(detsigma));
sic(i)=(-2*lik/T)+(NN*log(T))/T;
aic(i)=(-2*lik/T)+(2*NN)/T;
end
% sic
temp=1:maxlag;
if choice
 [trash,I]=min(aic);
else
   
[trash,I]=min(sic);
end

L=temp(I);