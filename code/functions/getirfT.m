function [ ir,hir,vir,dir ] = getirfT( L,LH,LV,N,horizon1,Fmat,Qmat,varcoef,iamat,A0,REPSx,Y,hlast,FF,scale,pos)
ir=0;
hir=0;
vir=0;
dir=0;
T=rows(Y);
hh=L+1:12:T;
parfor j=1:cols(hh)
   


Y0=Y(hh(j)-L+1:hh(j),:);
H0=hlast(hh(j)-LV+1:hh(j),1:N);
 
 [ir1,hir1,vir1,dir1] =...
    getlevelirfFMC( L,LH,LV,N,horizon1,Fmat,Qmat,varcoef,iamat,A0,REPSx,Y0,H0,FF,scale,pos );

ir=ir+ir1;
hir=hir+hir1;
vir=vir+vir1;
dir=dir+dir1;
end


hhx=cols(hh);
ir=ir/hhx;
hir=hir/hhx;
vir=vir/hhx;
dir=dir/hhx;