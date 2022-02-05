function [ out] = getvolirfFx( L,LH,LV,N,horizon,shock,Fmat,Qmat,varcoef )
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
        hhat(m,:)=[hx 0]*Fmat+vhat(m,:)*chol(Qmat(1:N,1:N));
        xhat=[];
        for i=1:L
            xhat=[xhat yhat(m-i,:)];
        end
        
        for i=1:LH
            xhat=[xhat hhat(m-i,:)];
        end
        xhat=[xhat 0];
        yhat(m,:)=xhat*varcoef;
    end
    
    out=[yhat(LL+1:end,:) hhat(LL+1:end,:)];
end

