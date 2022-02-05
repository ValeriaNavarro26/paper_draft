function [F1,F2,F3,F4,WW]=getspectrum(hhatx,L)
L=getlags(hhatx,L,2);
yy=hhatx;
N=cols(yy);
xx=[];
for j=1:L
xx=[xx lag0(hhatx,j)];
end
xx=[xx ones(rows(xx),1)];
yy=yy(L+1:end,:);
xx=xx(L+1:end,:);
b=xx\yy;
e=yy-xx*b;
sigx=(e'*e)/(rows(e)-cols(xx));
[F1,F2,F3,F4,WW]=spectrum(b,sigx,N,L);