function [ F11,F22,F33,F44,ww ] = ...
    spectrum( varcoef,varsigma,n,l )
TT=100;
TT2=ceil(TT/2);
%initialise frequency
w = zeros(TT,1);
ww=zeros(TT,1);
  im = sqrt(-1);
for j=1:TT
 w(j,1) = exp(-2.0*im*pi*(j-1)/TT); 
ww(j,1)=2*pi*(j-1)/TT;
end

F11=zeros(TT2,1);
F22=F11;
F33=F11;
F44=F11;

j=1;
    [FF,mu]=comp(vec(varcoef),n,l,1);
    
    OMEGA=zeros(rows(FF),rows(FF));
    OMEGA(1:n,1:n)=varsigma;
    for i=1:TT2
        tr=invpd(eye(cols(FF))-FF*w(i,1));
        g=(tr*OMEGA*tr')./(2*pi);
        
        %spectrum
        F11(i,j)=real(g(1,1));
        F22(i,j)=real(g(2,2));
        F33(i,j)=real(g(3,3));
        F44(i,j)=real(g(4,4));
      
        
    end

        