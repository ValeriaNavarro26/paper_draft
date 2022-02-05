function [ irx1,irx2,irx3,irx4,...
    virx1,virx2,virx3,virx4,...
    dirx1,dirx2,dirx3,dirx4,...
    irx1v,irx2v,irx3v,irx4v,...
    virx1v,virx2v,virx3v,virx4v,...
    dirx1v,dirx2v,dirx3v,dirx4v,...
    f1,f2,f3,f4,f5,f6,f7,f8,...
    fv1,fv2,fv3,fv4,fv5,fv6,fv7,fv8,...
    fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8] = getirfT4_0( L,LH,LV,N,horizon1,Fmat,Qmat,varcoef,iamat,REPSx,Y,hlast,FF)
irx1=0;
irx2=0;
irx3=0;
irx4=0;
irx5=0;

virx1=0;
virx2=0;virx3=0;virx4=0;virx5=0;
dirx1=0;dirx2=0;dirx3=0;dirx4=0;dirx5=0;

irx1v=0;
irx2v=0;
irx3v=0;
irx4v=0;
irx5v=0;

virx1v=0;
virx2v=0;virx3v=0;virx4v=0;virx5v=0;
dirx1v=0;dirx2v=0;dirx3v=0;dirx4v=0;dirx5v=0;
T=rows(Y);
hh=L+1:12:T;
parfor j=1:cols(hh)
%sigmax=iamat*diag(hlast(hh(j),1:N))*iamat';
A0=(iamat*diag(sqrt(hlast(hh(j),1:N))))';%chol(sigmax);


Y0=Y(hh(j)-L+1:hh(j),:);
H0=hlast(hh(j):hh(j),1:N);
 
 [ir1,ir2,ir3,ir4,hir1,hir2,hir3,hir4,...
    vir1,vir2,vir3,vir4,dir1,dir2,dir3,dir4] = ...
    getlevelirfFMCnew4_0( L,LH,LV,N,horizon1,Fmat,Qmat,varcoef,iamat,A0,REPSx,Y0,H0,FF);

irx1=irx1+ir1;irx2=irx2+ir2;irx3=irx3+ir3;irx4=irx4+ir4;

virx1=virx1+vir1;virx2=virx2+vir2;virx3=virx3+vir3;virx4=virx4+vir4;
dirx1=dirx1+dir1;dirx2=dirx2+dir2;dirx3=dirx3+dir3;dirx4=dirx4+dir4;

%volatility shocks


[ir1v,ir2v,ir3v,ir4v,hir1v,hir2v,hir3v,hir4v,...
    vir1v,vir2v,vir3v,vir4v,dir1v,dir2v,dir3v,dir4v] = ...
    getvolirfFMCnew4_0( L,LH,LV,N,horizon1,Fmat,Qmat,varcoef,iamat,A0,REPSx,Y0,H0,FF);

irx1v=irx1v+ir1v;irx2v=irx2v+ir2v;irx3v=irx3v+ir3v;irx4v=irx4v+ir4v;

virx1v=virx1v+vir1v;virx2v=virx2v+vir2v;virx3v=virx3v+vir3v;virx4v=virx4v+vir4v;
dirx1v=dirx1v+dir1v;dirx2v=dirx2v+dir2v;dirx3v=dirx3v+dir3v;dirx4v=dirx4v+dir4v;












end


hhx=cols(hh);
irx1=irx1/hhx;irx2=irx2/hhx;irx3=irx3/hhx;irx4=irx4/hhx;
virx1=virx1/hhx;virx2=virx2/hhx;virx3=virx3/hhx;virx4=virx4/hhx;
dirx1=dirx1/hhx;dirx2=dirx2/hhx;dirx3=dirx3/hhx;dirx4=dirx4/hhx;


irx1v=irx1v/hhx;irx2v=irx2v/hhx;irx3v=irx3v/hhx;irx4v=irx4v/hhx;
virx1v=virx1v/hhx;virx2v=virx2v/hhx;virx3v=virx3v/hhx;virx4v=virx4v/hhx;
dirx1v=dirx1v/hhx;dirx2v=dirx2v/hhx;dirx3v=dirx3v/hhx;dirx4v=dirx4v/hhx;





%compute FEVD
[ f1,f2,f3,f4,f5,f6,f7,f8] = getfvd4( irx1,irx2,irx3,irx4,irx1v,irx2v,irx3v,irx4v );
[ fv1,fv2,fv3,fv4,fv5,fv6,fv7,fv8 ] = getfvd4( virx1,virx2,virx3,virx4,virx1v,virx2v,virx3v,virx4v );
[ fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8 ] = getfvd4( dirx1,dirx2,dirx3,dirx4,dirx1v,dirx2v,dirx3v,dirx4v );







