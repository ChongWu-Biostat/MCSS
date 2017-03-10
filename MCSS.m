function beta_new = MCSS(A,lambda,tau1,tau2,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
% Input: 
%		A--n*p mutation matrix
%		lambda--a tuning parameter controlling the sparseness of the solutions
%		tau1--a tuning parameter in truncated lasso penalty
%		tau2--a tuning parameter relates to lambda
%		alpha--a tuning parameter relates to ridge penalty
%				The default value is 0.001.
%		beta_0--starting value of beta.
%				If we set beta_0=-1, a good starting value will be calculated.
%		method--method for step2 in Algorithm 1. 1 (default) is sub-gradient and  2  stands for CVX.
%		weight--weight of exclusivity cost. The default is 1.
%		max_iter_num--maximum number of iterations for difference of convex step.
%				The default value is 20.		
%		Err--termination condition for difference of convex step.
%				The default value is 0.001.
%		max_iter_num_s--maximum number of iterations for sub-gradient algorithm.
%				The default value is 20.
%		Err_s--termination condition for sub-gradient algorithm.
%				The default value is 0.001.
%
% Output: 
%		beta_new--a resulting estimate of beta.
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



[n,p]=size(A);Ap=sum(A,1);
beta_old=beta_0;    
if beta_0==-1
    beta_old=MCSS_0(A,weight,lambda,tau1,tau2);
end
err=2*Err+0.1;t=0;

while (t<=max_iter_num && err>Err)
    if method==1       % Subgradient
        beta_new_old=beta_old;
        err_s= 2*Err_s+0.1; tt=0;
        z_old2=Ap'.*(beta_old<=tau1)/tau1 ...
          +n*lambda*(beta_old<=tau2)/tau2-(1+weight)*Ap'/tau1;
        while (tt<=max_iter_num_s && err_s>Err_s)
            t_tt = 1/(2*sqrt(n*p)*(1+tt));
            g_old2 = z_old2+(1+weight)*A'*(A*beta_new_old/tau1>1)/(tau1)+2*alpha*beta_new_old;
            beta_new_new = beta_new_old-t_tt*g_old2;
            beta_new_new = beta_new_new.*(beta_new_new>=0);
            err_s = max(abs(beta_new_new-beta_new_old));beta_new_old=beta_new_new;tt=tt+1;
        end
        beta_new = beta_new_new;
    end
  
    if method==2       %CVX alpha>0
        cvx_begin quiet
      	variable beta_new(p);
      	minimize( alpha*beta_new'*beta_new+(1+weight)*sum(pos(A*beta_new/tau1-1))+beta_new'*(Ap'.*(beta_old<=tau1)/tau1+n*lambda*(beta_old<=tau2)/tau2-2*Ap'/tau1) );
     	subject to
         	tau1*ones(p,1)>=beta_new>=zeros(p,1); 
      cvx_end
    end

    err=sum(abs(beta_new-beta_old));
    beta_old=beta_new;
    t=t+1;
end
beta_new=beta_new.*(beta_new>1*1e-3);



