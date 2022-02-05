function [ out] = getlevelirfFx( L,LH,LV,N,horizon,shock,Fmat,Qmat,varcoef,A0 )

    LL=max(L,LH);
    hhat=zeros(horizon+LL,N);
    vhat=zeros(horizon+LL,N);
    vhat(LL+1,:)=shock;
    yhat=zeros(horizon+LL,N);
    for m=LL+1:horizon+LL
        hx=hhat(m-1,:);
        for i=1:LV
            hx=[hx yhat(m-i,:)];
        end
        hhat(m,:)=[hx 0]*Fmat;
        xhat=[];
        for i=1:L
            xhat=[xhat yhat(m-i,:)];
        end
        
        for i=1:LH
            xhat=[xhat hhat(m-i,:)];
        end
        xhat=[xhat 0];
        yhat(m,:)=xhat*varcoef+vhat(m,:)*A0;
    end
    
    out=[yhat(LL+1:end,:) hhat(LL+1:end,:)];
end

