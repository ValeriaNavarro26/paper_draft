function [ irx1,irx2,irx3,irx4,...
    virx1,virx2,virx3,virx4,...
    dirx1,dirx2,dirx3,dirx4,...
    irx1v,irx2v,irx3v,irx4v,...
    virx1v,virx2v,virx3v,virx4v,...
    dirx1v,dirx2v,dirx3v,dirx4v,...
    f1,f2,f3,f4,f5,f6,f7,f8,...
    fv1,fv2,fv3,fv4,fv5,fv6,fv7,fv8,...
    fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8] = ...
    getirfT4_tvp( L,LH,LV,N,horizon1,beta2v,Qmat,beta2,...
    iamat,REPSx,Y,hlast,T,Q,Qv,EX)
TT=floor(T/12);
irx1=zeros(TT,horizon1,N);
irx2=zeros(TT,horizon1,N);
irx3=zeros(TT,horizon1,N);
irx4=zeros(TT,horizon1,N);

virx1=zeros(TT,horizon1,N);
virx2=zeros(TT,horizon1,N);virx3=zeros(TT,horizon1,N);virx4=zeros(TT,horizon1,N);
dirx1=zeros(TT,horizon1,1);dirx2=zeros(TT,horizon1,1);
dirx3=zeros(TT,horizon1,1);dirx4=zeros(TT,horizon1,1);

irx1v=zeros(TT,horizon1,N);
irx2v=zeros(TT,horizon1,N);
irx3v=zeros(TT,horizon1,N);
irx4v=zeros(TT,horizon1,N);

virx1v=zeros(TT,horizon1,N);
virx2v=zeros(TT,horizon1,N);virx3v=zeros(TT,horizon1,N);virx4v=zeros(TT,horizon1,N);
dirx1v=zeros(TT,horizon1,1);dirx2v=zeros(TT,horizon1,1);
dirx3v=zeros(TT,horizon1,1);dirx4v=zeros(TT,horizon1,1);



f1=zeros(TT,horizon1,N);
f2=zeros(TT,horizon1,N);
f3=zeros(TT,horizon1,N);
f4=zeros(TT,horizon1,N);
f5=zeros(TT,horizon1,N);
f6=zeros(TT,horizon1,N);
f7=zeros(TT,horizon1,N);
f8=zeros(TT,horizon1,N);

fv1=zeros(TT,horizon1,N);
fv2=zeros(TT,horizon1,N);
fv3=zeros(TT,horizon1,N);
fv4=zeros(TT,horizon1,N);
fv5=zeros(TT,horizon1,N);
fv6=zeros(TT,horizon1,N);
fv7=zeros(TT,horizon1,N);
fv8=zeros(TT,horizon1,N);

fd1=zeros(TT,horizon1,1);
fd2=zeros(TT,horizon1,1);
fd3=zeros(TT,horizon1,1);
fd4=zeros(TT,horizon1,1);
fd5=zeros(TT,horizon1,1);
fd6=zeros(TT,horizon1,1);
fd7=zeros(TT,horizon1,1);
fd8=zeros(TT,horizon1,1);


T=rows(Y);
hh=L+1:12:T;
parfor j=1:cols(hh)

beta20=beta2(hh(j),:);
beta2v0=beta2v(hh(j),:);

A0=iamat';


Y0=Y(hh(j)-L+1:hh(j),:);
H0=hlast(hh(j):hh(j),1:N);
  
 [ir1,ir2,ir3,ir4,hir1,hir2,hir3,hir4,...
    vir1,vir2,vir3,vir4,dir1,dir2,dir3,dir4] = ...
    getlevelirfFMCnew4_tvp( L,LH,LV,N,horizon1,beta2v0,Qv,Qmat,...
    beta20,Q,iamat,A0,REPSx,Y0,H0,EX);





irx1(j,:,:)=ir1;irx2(j,:,:)=ir2;
irx3(j,:,:)=ir3;irx4(j,:,:)=ir4;


virx1(j,:,:)=vir1;virx2(j,:,:)=vir2;
virx3(j,:,:)=vir3;virx4(j,:,:)=vir4;

dirx1(j,:,:)=dir1;dirx2(j,:,:)=dir2;
dirx3(j,:,:)=dir3;dirx4(j,:,:)=dir4;












end














