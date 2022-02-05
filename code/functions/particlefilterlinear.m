function [ lik,states ] = ...
    particlefilterlinear( y,x,varcoef,iamat,F,cQ,FIT,...
    npart,T,N,L,EX,B00,P00)
%matrices of the state space
NS=cols(P00);


H=eye(N,NS);
lik=0;
states=zeros(T,NS);
%initial state
b0=repmat(B00,npart,1)+[randn(npart,N) zeros(npart,NS-N)]*chol(P00);
for i=1:T
    MU=zeros(1,NS);
    MU(1:N)=FIT(i,:);
    dens=zeros(npart,1);
	bnew=zeros(npart,NS);
    part=zeros(npart,NS);
    part(:,1:N)=randn(npart,N);
    for j=1:npart
        xtemp=zeros(1,N*L+EX);
        xtemp(:,1:N*L)=x(i,1:N*L);
        xtemp(:,end)=1;
        xtemp(:,(N*L)+1:end-1)=b0(j,N+1:end);
        %draw states
        htemp=MU+b0(j,:)*F'+part(j,:)*cQ;
        bnew(j,:)=htemp;
       
        
        %likelihood
        res=(y(i,:)-xtemp*varcoef);
        sigma=iamat*diag(exp(htemp(:,1:N)))*iamat';
        isigma=invpd(sigma);
        dsigma=det(sigma);
        
            resx=res*isigma*res';
            densx=(1/sqrt(dsigma))*exp(-0.5*resx);
       
        
        
	tempdens=densx;
    if isnan(tempdens) || isinf(tempdens) || ~isreal(tempdens)
        dens(j)=exp(-100);
    else
        dens(j)=tempdens;
    end
    end
    sdens=sum(dens);
   prob=dens./sdens;
   index1=discretesample(prob,npart);
   
   b0=bnew(index1,:);
   liki=sdens/npart;
lik=lik+log(liki);

states(i,:)=sum(b0.*repmat(prob,1,NS))';
end

