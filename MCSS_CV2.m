function  [beta_hat,M_list,CF_list] = MCSS_CV2(A,Lambda,Tau2, tau1,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s,M)
% Input:
%		A--n*p mutation matrix
%		Lambda--a set of tuning parameters for lambda.
%		Tau2--a set of tuning parameters for tau2.
%		tau1--a tuning parameter
%		alpha--a tuning parameter with ridge penalty
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
method = 1;
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

[n,p]=size(A);
TT = size(beta_0,2);

M_list=zeros(TT,p);
CF_list=zeros(TT,6);
beta_hat=zeros(p,1);
CF=Inf;
Ap=sum(A,1);
for tt=1:TT

	beta_0_tmp = beta_0(:,tt);
    [beta_old,~] = MCSS_CV1(A,Lambda,Tau2, tau1,alpha, beta_0_tmp, method,weight,max_iter_num, Err,max_iter_num_s,Err_s,M);
    M_tt=find(beta_old>0)';
    M_list(tt,1:length(M_tt))=M_tt;
    CF_new=-(1+weight)*sum(A*(beta_old>0)>0) + sum(A,1)*(beta_old>0);
    CF_list(tt,1)=CF_new;
    CF_list(tt,2)=-sum(A*(beta_old>0)>0)+sum(A,1)*(beta_old>0);
    CF_list(tt,3)=-weight*sum(A*(beta_old>0)>0);
    CF_list(tt,4)=sum(A*(beta_old>0)==1);
    CF_list(tt,5)=sum(A*(beta_old>0)==2);
    CF_list(tt,6)=sum(Ap(M_tt)<0.05*n);
    
    %[tt,CF_list(tt,1:5)/n]
    
    if CF>CF_new
        CF=CF_new;
        beta_hat=beta_old;
    end
end
CF_list=CF_list/n;


