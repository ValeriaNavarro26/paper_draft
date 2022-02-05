function [ir1,hir1,vir1,dir1] = getlevelirfFMC_3( L,LH,LV,N,horizon,Fmat,Qmat,varcoef,iamat,A0,REPS,Y0,hlast0,FF,scale,pos )
    %
    maxtrys=100;
    sigma0=A0'*A0;
     cQ=chol(Qmat);

    LL=max(L,LH);
    yy=0;
    yys=0;
    hh=0;
    hhs=0;
    vv=0;
    vvs=0;
dd=0;
dds=0;
    for R=1:REPS

    hhat=zeros(horizon+LL,N);
    hhat(L-LH+1:L,:)=log(hlast0);
    yhat=zeros(horizon+LL,N);
    yhat(1:L,:)=Y0;
    yhats=yhat;
    hhats=hhat;
    vol=zeros(horizon+LL,N);
    vols=vol;
detx=zeros(horizon+LL,1);
detxs=detx;
    for m=LL+1:horizon+LL
        hhatx=[];
        hhatxs=[];
        for i=1:LH
            hhatx=[hhatx hhat(m-i,:)];
            hhatxs=[hhatxs hhats(m-i,:)];
        end
   
        for i=1:LV
          hhatx=[hhatx yhat(m-i,:)];  
          hhatxs=[hhatxs yhats(m-i,:)]; 
        end
        
        if m==L+1
        hhat(m,:)=[hhatx 1]*Fmat;
        hhats(m,:)=[hhatxs 1]*Fmat;
        else
        chck=-1;
        trys=1;
        while chck<0 && trys<maxtrys;
        hhat(m,:)=[hhatx 1]*Fmat+randn(1,N)*cQ;
        hhats(m,:)=[hhatxs 1]*Fmat+randn(1,N)*cQ;
       
        ee1=exp(hhat(m,:))>2000;
        ee2=exp(hhats(m,:))>2000;
        ee=sum(ee1)+sum(ee2);
       
        if ee==0
            chck=1;
        else
            trys=trys+1;
        end
        if chck<0
            hhat(m,:)=log(hlast0(end,:));
            hhats(m,:)=log(hlast0(end,:));
        end
        end
        end
    % 
     %build covariance matrix
    sigma=iamat*diag(exp(hhat(m,:)))*iamat';
    sigmas=iamat*diag(exp(hhats(m,:)))*iamat';
    ee3=rcond(sigma);
    if ee3<1e-14 || isnan(ee3)
        sigma=sigma0;
    end
    ee4=rcond(sigmas);
    if ee4<1e-14 || isnan(ee4)
        sigmas=sigma0;
    end
    
    csigma=cholx(sigma);
    csigmas=cholx(sigmas);
    
   %simulate data wit & without a shock
      
        xhat=[];
        xhats=[];
        for i=1:L
            xhat=[xhat yhat(m-i,:)];
            xhats=[xhats yhats(m-i,:)];
        end
        
        for i=1:LH
            xhat=[xhat hhat(m-i,:)];
            xhats=[xhats hhats(m-i,:)];
        end
        xhat=[xhat 1];
        xhats=[xhats 1];
        
        
        if m==LL+1
shock=zeros(1,N);
shock(pos)=scale;
        yhats(m,:)=xhats*varcoef+shock*A0; 
        yhat(m,:)=xhat*varcoef;
        else
        yhats(m,:)=xhats*varcoef+randn(1,N)*csigmas; 
        yhat(m,:)=xhat*varcoef+randn(1,N)*csigma;
        end
      %unconditional volatility
    OMEGA=zeros(rows(FF),rows(FF));
    OMEGA(1:N,1:N)=sigma;
    uvar=doublej(FF,OMEGA);
    vol(m,:)=log(diag(uvar(1:N,1:N))');
    detx(m,:)=log(det(uvar(1:N,1:N)));
    OMEGAs=zeros(rows(FF),rows(FF));
    OMEGAs(1:N,1:N)=sigmas;
    uvars=doublej(FF,OMEGAs);
    vols(m,:)=log(diag(uvars(1:N,1:N))');
    detxs(m,:)=log(det(uvars(1:N,1:N)));    
        
    end
    
    yy=yy+yhat;
yys=yys+yhats;
hh=hh+hhat;
hhs=hhs+hhats;
vv=vv+vol;
vvs=vvs+vols;
dd=dd+detx;
dds=dds+detxs;
    end
 yy=yy/REPS;
yys=yys/REPS;
hh=hh/REPS;
hhs=hhs/REPS;
vv=vv/REPS;
vvs=vvs/REPS;
dd=dd/REPS;
dds=dds/REPS;

ir1=yys-yy;
hir1=hhs-hh;
vir1=vvs-vv;
dir1=dds-dd;

ir1=ir1(L+1:end,:);
hir1=hir1(L+1:end,:);
vir1=vir1(L+1:end,:);
dir1=dir1(L+1:end,:);