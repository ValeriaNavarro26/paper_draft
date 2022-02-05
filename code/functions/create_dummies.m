function [y,x] = create_dummies(lamda,tau,delta,epsilon,p,mu,sigma,n,ph,epsilonH)
% Creates matrices of dummy observations [...];
%lamda tightness parameter
%tau  prior on sum of coefficients
%delta prior mean for VAR coefficients
% epsilon tigtness of the prior around constant
% mu sample mean of the data
% sigma AR residual variances for the data




% Initialise output (necessary for final concatenation to work when tau=0):
x = [];
y = [];
yd1 = [];
yd2 = [];
xd1 = [];
xd2 = [];

%% Get dummy matrices in equation (5) of Banbura et al. 2007:

	yd1=[diag(sigma.*delta)./lamda;
         zeros(n*(p-1),n)];
     for i=1:ph*n
         yd1=[yd1;zeros(1,n)];
     end
     
         
        yd1=[yd1; zeros(1,n)];
     
	jp=diag(1:p);
    
	xd1=[kron(jp,diag(sigma)./lamda)];
  
  if epsilonH>0
    for i=1:ph*n
         xd1=[[xd1 zeros(rows(xd1),1)]; [zeros(1,cols(xd1)) 1./sqrt(epsilonH)]];
     end
  end
         xd1=[[xd1 zeros(rows(xd1),1)]; [zeros(1,cols(xd1)) epsilon]];
   


%% Get additional dummy matrices - see equation (9) of Banbura et al. 2007:
if tau>0
    
	yd2=diag(delta.*mu)./tau;
	xd2=[kron(ones(1,p),yd2) zeros(n,n*ph) zeros(n,1)];
    
end
     
%% 
y=[yd1;yd2];
% xd1
% size(xd1)
% size(xd2)
% size(kron(ones(1,p),yd2))
x=[xd1;xd2];
 
         
 
 