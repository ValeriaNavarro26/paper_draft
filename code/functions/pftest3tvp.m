
function [beta,BB,W]=pftest3tvp(beta2v,cQ,iQ,iamat,xfix,y,x,N,NS,B00,beta2,Nx,LV,EX,L,X0)
%




T=size(y,1);


beta=zeros(T,NS,N); %will hold filtered state variable
beta(1,:,:)=repmat(B00',1,N); %fixed initial condition
BB=zeros(N,T); %ancestor indices 
W=zeros(N,T); %weights
beta(1,:,N)=xfix(1,:);
for t=1:T
    %matrices of the state space
    betaf1=reshape(beta2v(t,:),Nx*LV+(Nx+1),Nx);
F=zeros(NS,NS);
F(1:Nx,1:Nx)=betaf1(1:Nx,:)';F(Nx+1:NS,1:NS-Nx)=eye(NS-Nx);
fit=X0(t,Nx+1:end)*betaf1(Nx+1:end,:);
varcoef=reshape(beta2(t,:),Nx*L+EX,Nx);
    
if (t~=1)
    index=discretesample(W(:,t-1),N);
    index=index(randperm(N));
 
    [betahat,m]=xpred(beta(t-1,:,:),F,fit,xfix(t,:),iQ);
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
R=iamat*diag(exp(beta1))*iamat';
iR=invpd(R);
detR=det(R);
res=y(t,:)-x(t,:)*varcoef;
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
for j=1:N
 fitx= squeeze(beta(:,:,j))*F'+fit;
 res=fitx-x;
 densx=exp(-0.5*res*iQ*res');
 %squeeze(betahat(:,:,index))'+randn(N,NS)'*cQ
out(:,:,j)=fitx;
dens(j,:)=densx;
end

end