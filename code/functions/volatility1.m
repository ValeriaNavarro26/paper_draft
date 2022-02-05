function [pred] = volatility1( ss,iamat,hh,n,l,ex )

TT=size(hh,1);
pred=zeros(TT,n);
[FF,mu]=comp(ss,n,l,ex);
    iA=iamat;
for j=1:TT
    
    
    covmat=iA*diag(hh(j,:))*iA';
    OMEGA=zeros(rows(FF),rows(FF));
    OMEGA(1:n,1:n)=covmat;
    %unconditional variance
    uvar=invpd(eye(cols(FF),cols(FF))-kron(FF,FF))*vec(OMEGA);
    uvar=reshape(uvar,cols(FF),cols(FF));
    
    pred(j,:)=diag(uvar(1:n,1:n))';
    
    
end