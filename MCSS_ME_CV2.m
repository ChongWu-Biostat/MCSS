function  [beta_hat,M_list,CF_list] = MCSS_ME_CV2(A,G,gamma,Lambda,Tau2, tau1,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s,M)
% Input:
%		A--n*p mutation matrix
%		G--n*p gene expression matrix
%		gamma--a tuning parameter controlling the contributions between mutation data and gene expression data
%		Lambda--a set of tuning parameters for lambda.
%		Tau2--a set of tuning parameters for tau2.
%		tau1--tuning parameters
%		alpha--tuning parameters with ridge penalty
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
%		M_list--resulting list of gene sets corresponding to Lambda
%		CF--cost values corresponding to Lambda
%
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0


if ~exist('gamma', 'var')
    Err_s = 0.5;
end

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
CF_list=zeros(TT,7);
beta_hat=zeros(p,1);
CF=Inf;vn=20;

Ap=sum(A,1);[~,b4]=sort(Ap);b4=b4(max(1,(p-vn+1)):p);


for tt=1:TT
    beta_0_tmp = beta_0(:,tt);
	[beta_old,~] = MCSS_ME_CV1(A,G, gamma, Lambda,Tau2, tau1,alpha, beta_0_tmp, method,weight,max_iter_num, Err,max_iter_num_s,Err_s,M);
    
    M_tt=find(beta_old>0)';
    M_list(tt,1:length(M_tt))=M_tt;Cor=abs(corr(G));Cor=Cor-eye(p);
    CF_new=-(1+weight)*sum(A*(beta_old>0)>0)/n+sum(A,1)*(beta_old>0)/n-gamma*(beta_old>0)'*Cor*(beta_old>0)/(sum((beta_old>0))^2+1e-5);
    CF_list(tt,1)=CF_new;
    CF_list(tt,2)=-sum(A*(beta_old>0)>0)/n+sum(A,1)*(beta_old>0)/n;
    CF_list(tt,3)=-weight*sum(A*(beta_old>0)>0)/n;
    CF_list(tt,4)=sum(A*(beta_old>0)==1)/n;
    CF_list(tt,5)=sum(A*(beta_old>0)==2)/n;
    CF_list(tt,6)=sum(Ap(M_tt)<0.05*n)/n;
    CF_list(tt,7)=-gamma*(beta_old>0)'*Cor*(beta_old>0)/(sum((beta_old>0))^2+1e-5);
    %[tt,CF_list(tt,1:7)]
    if CF>CF_new
        CF=CF_new;
        beta_hat = beta_old;
    end
end

