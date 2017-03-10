function beta_0 = MCSS_0(A,weight,Lambda,tau1,tau2)
% Input: A--n*p mutation matrix
%        weigth--weigth of exclusivity cost
%        Lambda--tuning set of lambda
%        tau1,tau2--tuning parameters
% Output: beta_0--an initial estimate of beta for MCSS.
[n,p]=size(A);
Ap=sum(A,1);
CF=Inf;
for  i=1:length(Lambda)
    lambda=Lambda(i);
    f=[n*lambda/tau2+Ap'/tau1; -(1+weight)*ones(n,1)];
    a=[-A/tau1, eye(n)];
    b=zeros(n,1);
    lb=zeros(p+n,1);
    ub=[tau1*ones(p,1); ones(n,1)];
    beta_old=linprog(f,a,b,[],[],lb,ub); 
    beta_old=beta_old(1:p,1);
    CF_new=-(1+weight)*sum(A*(beta_old>0)>0)+sum(A,1)*(beta_old>0);
    if CF>CF_new
        CF=CF_new;
        beta_0=beta_old;
    end
end
