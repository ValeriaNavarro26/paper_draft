function [ outL,outh,tmp,tmpv ] = getirfFx(L,LH,N,LV,horizon1,Fmat,Qmat,varcoef,A0  )
mse=0;
for i=1:N
    shock=zeros(1,N);
    shock(i)=1;
      tmp= getlevelirfFx( L,LH,LV,N,horizon1,shock,Fmat,Qmat,varcoef,A0 );
       tmpv= getvolirfFx( L,LH,LV,N,horizon1,shock,Fmat,Qmat,varcoef );
       outL(i,:,:)=tmp;
       outh(i,:,:)=tmpv;
       mse=mse+(cumsum(tmp.^2)+cumsum(tmpv.^2));
       
      
end



