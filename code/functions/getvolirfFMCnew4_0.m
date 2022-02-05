function [ir1,ir2,ir3,ir4,hir1,hir2,hir3,hir4,...
    vir1,vir2,vir3,vir4,dir1,dir2,dir3,dir4] = ...
    getvolirfFMCnew4_0( L,LH,LV,N,...
    horizon,Fmat,Qmat,varcoef,iamat,A0,REPS,Y0,hlast0,FF )
    %
    maxtrys=100;
    sigma0=A0'*A0;
     cQ=chol(Qmat);

    LL=max(L,LH);
    yy=0;
    yys1=0;yys2=0;yys3=0;yys4=0;
    hh=0;
    hhs1=0;hhs2=0;hhs3=0;hhs4=0;
    vv=0;
    vvs1=0;vvs2=0;vvs3=0;vvs4=0;
dd=0;
dds1=0;dds2=0;dds3=0;dds4=0;
    for R=1:REPS

    hhat=zeros(horizon+LL,N);
    hhat(LL:LL,:)=log(hlast0);
    yhat=zeros(horizon+LL,N);
    yhat(1:L,:)=Y0;
    yhats1=yhat;
    hhats1=hhat;
    yhats2=yhat;
    hhats2=hhat;
    yhats3=yhat;
    hhats3=hhat;
    yhats4=yhat;
    hhats4=hhat;
 
    
    vol=zeros(horizon+LL,N);
    vols1=vol;
     vols2=vol;
      vols3=vol;
       vols4=vol;
       
detx=zeros(horizon+LL,1);
detxs1=detx;
detxs2=detx;
detxs3=detx;
detxs4=detx;

    for m=LL+1:horizon+LL
        hhatx=hhat(m-1,:);
        hhatxs1=hhats1(m-1,:);
        hhatxs2=hhats2(m-1,:);
        hhatxs3=hhats3(m-1,:);
        hhatxs4=hhats4(m-1,:);
    
        for i=1:LV
              hhatx=[hhatx yhat(m-i,:)];  
          hhatxs1=[hhatxs1 yhats1(m-i,:)]; 
          hhatxs2=[hhatxs2 yhats2(m-i,:)]; 
          hhatxs3=[hhatxs3 yhats3(m-i,:)]; 
          hhatxs4=[hhatxs4 yhats4(m-i,:)];
          
        end
        if m==L+1
        hhat(m,:)=[hhatx 1]*Fmat;
       shock=zeros(1,N);shock(1)=1; hhats1(m,:)=[hhatxs1 1]*Fmat+shock*cQ;
       shock=zeros(1,N);shock(2)=1;hhats2(m,:)=[hhatxs2 1]*Fmat+shock*cQ;
       shock=zeros(1,N);shock(3)=1;hhats3(m,:)=[hhatxs3 1]*Fmat+shock*cQ;
       shock=zeros(1,N);shock(4)=1;hhats4(m,:)=[hhatxs4 1]*Fmat+shock*cQ;
        
        else
        chck=-1;
        trys=1;
        while chck<0 && trys<maxtrys;
        hhat(m,:)=[hhatx 1]*Fmat+randn(1,N)*cQ;
        hhats1(m,:)=[hhatxs1 1]*Fmat+randn(1,N)*cQ;
        hhats2(m,:)=[hhatxs2 1]*Fmat+randn(1,N)*cQ;
        hhats3(m,:)=[hhatxs3 1]*Fmat+randn(1,N)*cQ;
        hhats4(m,:)=[hhatxs4 1]*Fmat+randn(1,N)*cQ;
        
       
        ee1=exp(hhat(m,:))>2000;
        ee2=exp(hhats1(m,:))>2000;
        ee3=exp(hhats2(m,:))>2000;
        ee4=exp(hhats3(m,:))>2000;
        ee5=exp(hhats4(m,:))>2000;
        
        ee=sum(ee1)+sum(ee2)+sum(ee3)+sum(ee4)+sum(ee5);
       
        if ee==0
            chck=1;
        else
            trys=trys+1;
        end
        if chck<0
            hhat(m,:)=log(hlast0(end,:));
            hhats1(m,:)=log(hlast0(end,:));
            hhats2(m,:)=log(hlast0(end,:));
            hhats3(m,:)=log(hlast0(end,:));
            hhats4(m,:)=log(hlast0(end,:));
            
        end
        end
        end
    % 
     %build covariance matrix
    sigma=iamat*diag(exp(hhat(m,:)))*iamat';
    sigmas1=iamat*diag(exp(hhats1(m,:)))*iamat';
    sigmas2=iamat*diag(exp(hhats2(m,:)))*iamat';
    sigmas3=iamat*diag(exp(hhats3(m,:)))*iamat';
    sigmas4=iamat*diag(exp(hhats4(m,:)))*iamat';
    
    ee7=rcond(sigma);
    if ee7<1e-14 || isnan(ee7)
        sigma=sigma0;
    end
    ee8=rcond(sigmas1);
    if ee8<1e-14 || isnan(ee8)
        sigmas1=sigma0;
    end
    
    ee9=rcond(sigmas2);
    if ee9<1e-14 || isnan(ee9)
        sigmas2=sigma0;
    end
    
    ee10=rcond(sigmas3);
    if ee10<1e-14 || isnan(ee10)
        sigmas3=sigma0;
    end
    
    ee11=rcond(sigmas4);
    if ee11<1e-14 || isnan(ee11)
        sigmas4=sigma0;
    end
    
    
    
    csigma=cholx(sigma);
    csigmas1=cholx(sigmas1);
    csigmas2=cholx(sigmas2);
    csigmas3=cholx(sigmas3);
    csigmas4=cholx(sigmas4);
    
   %simulate data wit & without a shock
      
        xhat=[];
        xhats1=[];xhats2=[];xhats3=[];xhats4=[];
        for i=1:L
            xhat=[xhat yhat(m-i,:)];
            xhats1=[xhats1 yhats1(m-i,:)];
            xhats2=[xhats2 yhats2(m-i,:)];
            xhats3=[xhats3 yhats3(m-i,:)];
            xhats4=[xhats4 yhats4(m-i,:)];
            
        end
        
        for i=1:LH
            xhat=[xhat hhat(m-i,:)];
            xhats1=[xhats1 hhats1(m-i,:)];
            xhats2=[xhats2 hhats2(m-i,:)];
            xhats3=[xhats3 hhats3(m-i,:)];
            xhats4=[xhats4 hhats4(m-i,:)];
            
        end
        xhat=[xhat 1];
        xhats1=[xhats1 1];
        xhats2=[xhats2 1];
        xhats3=[xhats3 1];
        xhats4=[xhats4 1];
        
        
       
        yhats1(m,:)=xhats1*varcoef+randn(1,N)*csigmas1; 
        yhats2(m,:)=xhats2*varcoef+randn(1,N)*csigmas2;
        yhats3(m,:)=xhats3*varcoef+randn(1,N)*csigmas3;
        yhats4(m,:)=xhats4*varcoef+randn(1,N)*csigmas4;
        
        yhat(m,:)=xhat*varcoef+randn(1,N)*csigma;
       
      %unconditional volatility
    OMEGA=zeros(rows(FF),rows(FF));
    OMEGA(1:N,1:N)=sigma;
    uvar=doublej(FF,OMEGA);
    vol(m,:)=log(diag(uvar(1:N,1:N))');
    detx(m,:)=log(det(uvar(1:N,1:N)));
    %
    OMEGAs1=zeros(rows(FF),rows(FF));
    OMEGAs1(1:N,1:N)=sigmas1;
    uvars1=doublej(FF,OMEGAs1);
    vols1(m,:)=log(diag(uvars1(1:N,1:N))');
    detxs1(m,:)=log(det(uvars1(1:N,1:N)));   
    %
      OMEGAs2=zeros(rows(FF),rows(FF));
    OMEGAs2(1:N,1:N)=sigmas2;
    uvars2=doublej(FF,OMEGAs2);
    vols2(m,:)=log(diag(uvars2(1:N,1:N))');
    detxs2(m,:)=log(det(uvars2(1:N,1:N))); 
    %
        OMEGAs3=zeros(rows(FF),rows(FF));
    OMEGAs3(1:N,1:N)=sigmas3;
    uvars3=doublej(FF,OMEGAs3);
    vols3(m,:)=log(diag(uvars3(1:N,1:N))');
    detxs3(m,:)=log(det(uvars3(1:N,1:N))); 
    %
        OMEGAs4=zeros(rows(FF),rows(FF));
    OMEGAs4(1:N,1:N)=sigmas4;
    uvars4=doublej(FF,OMEGAs4);
    vols4(m,:)=log(diag(uvars4(1:N,1:N))');
    detxs4(m,:)=log(det(uvars4(1:N,1:N))); 
    %
       
        
    end
    
    yy=yy+yhat;
yys1=yys1+yhats1;
yys2=yys2+yhats2;
yys3=yys3+yhats3;
yys4=yys4+yhats4;

hh=hh+hhat;
hhs1=hhs1+hhats1;
hhs2=hhs2+hhats2;
hhs3=hhs3+hhats3;
hhs4=hhs4+hhats4;
vv=vv+vol;
vvs1=vvs1+vols1;
vvs2=vvs2+vols2;
vvs3=vvs3+vols3;
vvs4=vvs4+vols4;

dd=dd+detx;
dds1=dds1+detxs1;
dds2=dds2+detxs2;
dds3=dds3+detxs3;
dds4=dds4+detxs4;

    end
 yy=yy/REPS;
yys1=yys1/REPS;yys2=yys2/REPS;yys3=yys3/REPS;yys4=yys4/REPS;
hh=hh/REPS;
hhs1=hhs1/REPS;hhs2=hhs2/REPS;hhs3=hhs3/REPS;hhs4=hhs4/REPS;
vv=vv/REPS;
vvs1=vvs1/REPS;vvs2=vvs2/REPS;vvs3=vvs3/REPS;vvs4=vvs4/REPS;
dd=dd/REPS;
dds1=dds1/REPS;dds2=dds2/REPS;dds3=dds3/REPS;dds4=dds4/REPS;

ir1=yys1-yy;ir2=yys2-yy;ir3=yys3-yy;ir4=yys4-yy;
hir1=hhs1-hh;hir2=hhs2-hh;hir3=hhs3-hh;hir4=hhs4-hh;
vir1=vvs1-vv;vir2=vvs2-vv;vir3=vvs3-vv;vir4=vvs4-vv;
dir1=dds1-dd;dir2=dds2-dd;dir3=dds3-dd;dir4=dds4-dd;

ir1=ir1(L+1:end,:);ir2=ir2(L+1:end,:);ir3=ir3(L+1:end,:);
ir4=ir4(L+1:end,:);
hir1=hir1(L+1:end,:);hir2=hir2(L+1:end,:);hir3=hir3(L+1:end,:);
hir4=hir4(L+1:end,:);
vir1=vir1(L+1:end,:);vir2=vir2(L+1:end,:);vir3=vir3(L+1:end,:);
vir4=vir4(L+1:end,:);
dir1=dir1(L+1:end,:);dir2=dir2(L+1:end,:);dir3=dir3(L+1:end,:);
dir4=dir4(L+1:end,:);