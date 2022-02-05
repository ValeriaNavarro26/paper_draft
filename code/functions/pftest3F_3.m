
function [beta,BB,W]=pftest3F_3(F,fit,cQ,iQ,iamat,xfix,y,x,N,NS,B00,B,Nx)
%




T=size(y,1);


beta=zeros(T,NS,N); %will hold filtered state variable
beta(1,:,:)=repmat(B00',1,N); %fixed initial condition
BB=zeros(N,T); %ancestor indices 
W=zeros(N,T); %weights
beta(1,:,N)=xfix(1,:);
for t=1:T
if (t~=1)
    index=discretesample(W(:,t-1),N);
    index=index(randperm(N));
 
    [betahat,m]=xpred(beta(t-1,:,:),F,fit(t,:),xfix(t,:),iQ);
    %draw particles
    for j=1:N
    beta(t,:,j)=squeeze(betahat(:,:,index(j)))+randn(1,NS)*cQ;
    end
    %conditioned particle
    beta(t,:,N)=xfix(t,:);
    %sample ancestors
    was=m.*W(:,t-1);
    probas=was./sum(was);
    index(N)=discretesample(probas,1);%find(rand(1) < cumsum(probas),1,'first');%discretesample(probas,1);

    BB(:,t)=index';
end

%compute importance weights
logweights=zeros(1,N);
for j=1:N
beta1=squeeze(beta(t,:,j));
R=iamat*diag(exp(beta1(1:Nx)))*iamat';
iR=invpd(R);
detR=det(R);
xvar=[x(t,1:cols(x)-1) beta1(Nx+1:end) x(t,end)];

res=y(t,:)-xvar*B;
res2=res*iR*res';
logweights(j)=log((1/sqrt(detR)))+(-0.5*res2);
end
const = max(logweights); % Subtract the maximum value for numerical stability
weights = exp(logweights-const);
W(:,t)=weights./sum(weights);


end


% Generate the trajectories from ancestor indices
index = BB(:,T);
for t = T-1:-1:1
    beta(t,:,:) = beta(t,:,index);
    index = BB(index,t);
end

end


function [out,dens]= xpred(beta,F,fit,x,iQ)
N=size(beta,3);
NS=size(beta,2);
out=zeros(1,NS,N);
dens=zeros(N,1);
mu=zeros(1,NS);
mu(:,1:cols(fit))=fit;
for j=1:N
 fitx= squeeze(beta(:,:,j))*F'+mu;
 res=fitx-x;
 densx=exp(-0.5*res*iQ*res');
 %squeeze(betahat(:,:,index))'+randn(N,NS)'*cQ
out(:,:,j)=fitx;
dens(j,:)=densx;
end

end