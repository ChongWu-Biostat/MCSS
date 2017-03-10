function [beta_t,CF] = MCSS_CV1(A,Lambda,Tau2, tau1,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s,M)
% Input:
%		A--n*p mutation matrix
%		Lambda--a set of tuning parameters for lambda.
%		Tau2--a set of tuning parameters for tau2.
%		tau1--a tuning parameter in truncated lasso penalty
%		alpha--a tuning parameter relates to ridge penalty
%				The default value is 0.001.
%		beta_0--starting value of beta.
%				If we set beta_0=-1, a good starting value will be calculated.
%		method--method for step2 in Algorithm 1. 2 (default) stands for CVX.
%		weight--weight of exclusivity cost. The default is 1.
%		max_iter_num--maximum number of iterations for difference of convex step.
%				The default value is 20.		
%		Err--termination condition for difference of convex step.
%				The default value is 0.001.
%		max_iter_num_s--maximum number of iterations for sub-gradient algorithm.
%				The default value is 20.
%		Err_s--termination condition for sub-gradient algorithm.
%				The default value is 0.001.
%		M--repeat CV for M times
%
% Output: 
%		beta_t--a resulting estimate of beta
%       CF--cost values corresponding Lambda
%
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0


if ~exist('Err_s', 'var')
    Err_s = 1*1e-3;
end


if ~exist('max_iter_num_s', 'var')
    max_iter_num_s = 20;
end


if ~exist('Err', 'var')
    Err = 1*1e-3;
end

if ~exist('max_iter_num', 'var')
    max_iter_num = 20;
end

if ~exist('weight', 'var')
    weight = 1;
end

if ~exist('method', 'var')
    method = 2;
end

if ~exist('beta_0', 'var')
    beta_0 = -1;
end

if ~exist('alpha', 'var')
    alpha=1*1e-3;
end

if ~exist('M', 'var')
    M=1;
end



[n,~]=size(A);
K=length(Lambda);
KK=length(Tau2);
if (K+KK)==2
    beta_t = MCSS(A,Lambda,tau1,Tau2,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s);
end
if (K+KK)~=2    
    n_tr=(n-mod(n,2))/2;
    CF=zeros(K,KK);
    for m=1:M
        s=RandStream('mcg16807','Seed',101);
        RandStream.setGlobalStream(s);
        perm=randperm(n);
        A_perm=A(perm,:); 
        A_tr=A_perm(1:n_tr,:);  
        A_tu=A_perm((n_tr+1):n,:);  
        for k=1:K
          lambda=Lambda(k);
          for kk=1:KK
            tau2=Tau2(kk)*tau1; 
            beta_k = MCSS(A_tr,lambda,tau1,tau2,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s);
            CF(k,kk)=CF(k,kk)+((1+weight)*sum(A_tu*(beta_k>0)>0)-sum(A_tu,1)*(beta_k>0));
          end
        end
    end
    CF=CF/M;
    CF_max=max(max(CF));
    i_CF_max=find(CF==CF_max);
    [Ib,Jb]=ind2sub(size(CF),i_CF_max(1));
    tau2_hat=tau1*Tau2(Jb);
    lambda_hat=Lambda(Ib);
    beta_t = MCSS(A,lambda_hat,tau1,tau2_hat,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s);
end
