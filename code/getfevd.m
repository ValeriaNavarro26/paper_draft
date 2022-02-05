clear
addpath('C:\Users\Valeria\Desktop\CICLO 2020-1\TESIS\Benchmarkmodel\Commodity\functions')
parpool(6)
load priors
horizon1=60;
scale=1;

 load(sfile,'amatsave','bsave','fsave','qsave','musave','hsave','naccept');
fsize=500;
REPSx=50;
irfmat1=zeros(fsize,horizon1,N*2+1);
vmat1=zeros(fsize,horizon1,N*2+1);
irfmat2=zeros(fsize,horizon1,N*2+1);
vmat2=zeros(fsize,horizon1,N*2+1);
irfmat3=zeros(fsize,horizon1,N*2+1);
vmat3=zeros(fsize,horizon1,N*2+1);
irfmat4=zeros(fsize,horizon1,N*2+1);
vmat4=zeros(fsize,horizon1,N*2+1);
irfmat5=zeros(fsize,horizon1,N*2+1);
vmat5=zeros(fsize,horizon1,N*2+1);
%
irfmat1v=zeros(fsize,horizon1,N*2+1);
vmat1v=zeros(fsize,horizon1,N*2+1);
irfmat2v=zeros(fsize,horizon1,N*2+1);
vmat2v=zeros(fsize,horizon1,N*2+1);
irfmat3v=zeros(fsize,horizon1,N*2+1);
vmat3v=zeros(fsize,horizon1,N*2+1);
irfmat4v=zeros(fsize,horizon1,N*2+1);
vmat4v=zeros(fsize,horizon1,N*2+1);
irfmat5v=zeros(fsize,horizon1,N*2+1);
vmat5v=zeros(fsize,horizon1,N*2+1);



for jgibbs=1:fsize
    iamat=squeeze(amatsave(jgibbs,:,:));
    beta2=bsave(jgibbs,:,:);
fbig=squeeze(fsave(jgibbs,:,:))';
Qbig=squeeze(qsave(jgibbs,:,:))';
mubig=squeeze(musave(jgibbs,:,:))';
hlast=squeeze(hsave(jgibbs,:,:)); 
betaf1=reshape(fbig,N+(N*LV+1),N);
Fmat=betaf1;

Qmat=reshape(Qbig,N,N);
mumat=zeros(N,1);
for j=1:N

mumat(j,1)=mubig(1,j);
end
FF=comp(beta2,N,L,EX);
varcoef=reshape(beta2,N*L+EX,N);



 jgibbs
 [ irx1,irx2,irx3,irx4,...
    virx1,virx2,virx3,virx4,...
    dirx1,dirx2,dirx3,dirx4,...
    irx1v,irx2v,irx3v,irx4v,...
    virx1v,virx2v,virx3v,virx4v,...
    dirx1v,dirx2v,dirx3v,dirx4v,...
    f1,f2,f3,f4,f5,f6,f7,f8,...
    fv1,fv2,fv3,fv4,fv5,fv6,fv7,fv8,...
    fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8] = getirfT4( L,LH,LV,N,horizon1,Fmat,Qmat,varcoef,iamat,REPSx,y,hlast(2:end,:),FF);
tmp1=[irx1 virx1 dirx1];
tmp1f=[f1 fv1 fd1];
tmp1v=[irx1v virx1v dirx1v];
tmp1fv=[f5 fv5 fd5];
irfmat1(jgibbs,:,:)=tmp1;
vmat1(jgibbs,:,:)=tmp1f;
irfmat1v(jgibbs,:,:)=tmp1v;
vmat1v(jgibbs,:,:)=tmp1fv;
%
tmp1=[irx2 virx2 dirx2];
tmp1f=[f2 fv2 fd2];
tmp1v=[irx2v virx2v dirx2v];
tmp1fv=[f6 fv6 fd6];
irfmat2(jgibbs,:,:)=tmp1;
vmat2(jgibbs,:,:)=tmp1f;
irfmat2v(jgibbs,:,:)=tmp1v;
vmat2v(jgibbs,:,:)=tmp1fv;
%
tmp1=[irx3 virx3 dirx3];
tmp1f=[f3 fv3 fd3];
tmp1v=[irx3v virx3v dirx3v];
tmp1fv=[f7 fv7 fd7];
irfmat3(jgibbs,:,:)=tmp1;
vmat3(jgibbs,:,:)=tmp1f;
irfmat3v(jgibbs,:,:)=tmp1v;
vmat3v(jgibbs,:,:)=tmp1fv;
%
tmp1=[irx4 virx4 dirx4];
tmp1f=[f4 fv4 fd4];
tmp1v=[irx4v virx4v dirx4v];
tmp1fv=[f8 fv8 fd8];
irfmat4(jgibbs,:,:)=tmp1;
vmat4(jgibbs,:,:)=tmp1f;
irfmat4v(jgibbs,:,:)=tmp1v;
vmat4v(jgibbs,:,:)=tmp1fv;
%
    
    
   
end
save('fevdtestnew');



VAL=prctile(vmat1v,[50 5 95 16 84 ],1);

MP_4_12=zeros(12,1)
for j=1:4
MP_4_12(j,1)=VAL(1,12,j)
MP_4_12(4+j,1)=VAL(1,24,j)
MP_4_12(8+j,1)=VAL(1,60,j)
end

PBI_4_12=zeros(12,1)
for j=1:4
PBI_4_12(j,1)=VAL(1,12,j)
PBI_4_12(4+j,1)=VAL(1,24,j)
PBI_4_12(8+j,1)=VAL(1,60,j)
end

IPC_4_12=zeros(12,1)
for j=1:4
IPC_4_12(j,1)=VAL(1,12,j)
IPC_4_12(4+j,1)=VAL(1,24,j)
IPC_4_12(8+j,1)=VAL(1,60,j)
end

TC_4_12=zeros(12,1)
for j=1:4
TC_4_12(j,1)=VAL(1,12,j)
TC_4_12(4+j,1)=VAL(1,24,j)
TC_4_12(8+j,1)=VAL(1,60,j)
end

VAL=prctile(vmat1v,[50 5 95 16 84 ],1);

MP_4_LB=zeros(12,1)
for j=1:4
    MP_4_LB(j,1)=VAL(5,12,j)
    MP_4_LB(4+j,1)=VAL(5,12,j)
    MP_4_LB(8+j,1)=VAL(5,12,j)
end

VAL=prctile(vmat2v,[50 5 95 16 84 ],1);

PBI_4_LB=zeros(12,1)
for j=1:4
    PBI_4_LB(j,1)=VAL(5,12,j)
    PBI_4_LB(4+j,1)=VAL(5,12,j)
    PBI_4_LB(8+j,1)=VAL(5,12,j)
end
VAL=prctile(vmat3v,[50 5 95 16 84 ],1);
IPC_4_LB=zeros(12,1)
for j=1:4
    IPC_4_LB(j,1)=VAL(5,12,j)
    IPC_4_LB(4+j,1)=VAL(5,12,j)
    IPC_4_LB(8+j,1)=VAL(5,12,j)
end


VAL=prctile(vmat4v,[50 5 95 16 84 ],1);
TC_4_LB=zeros(12,1)
for j=1:4
    TC_4_LB(j,1)=VAL(5,12,j)
    TC_4_LB(4+j,1)=VAL(5,12,j)
    TC_4_LB(8+j,1)=VAL(5,12,j)
end





