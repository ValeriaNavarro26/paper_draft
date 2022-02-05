function [beta,res,problem,trys]=getvarcoef2(y,x,yd,xd,Sbig,iamat,lamda,maxtrys,L,EX,iH0B0,iH0);

%normalise
[T,N]=size(y);
lamday=repmat(sqrt(lamda),1,N);
lamdax=repmat(sqrt(lamda),1,cols(x));
Y=y./lamday;
X=x./lamdax;
Y0=[Y;yd];
X0=[X;xd];
%conditional mean of the VAR coefficients
 bhat=vec(X0\Y0);  %ols on the appended data
 xx=(X0'*X0);
 sigma=iamat*diag(Sbig)*iamat';
 vstar=invpd(iH0+kron(invpd(sigma),xx));
 mstar=vstar*(iH0B0+kron(invpd(sigma),xx)*bhat);

trys=1;
problem=0;
chck=-1;
while chck<0 && trys<maxtrys
beta=mstar+(randn(1,N*(N*L+EX))*cholx(vstar))';
ee=stability(beta,N,L,EX);
if ee==0;
    chck=10;
else
    trys=trys+1;
end
end
if chck<0
    problem=1;
end
res=y-x*reshape(beta,N*L+EX,N);