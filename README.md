#MCSS  

Cancer is characterized by numerous somatic mutations, however, only a subset of mutations, called *driver* mutations, contributes to tumor growth and progression. 

Recent studies indicated that mutations arising in a driver pathway often covers most samples and for a single sample only a singe or few mutations appear because a single mutation is capable to perturb the whole pathway. 

Based on these findings, minimum cost subset selection (**MCSS**) is developed for de novo discovery of mutated driver pathway in cancer studies. 

## Features

* We provide a novel regularization based way (MCSS) for de novo discovery mutated driver pathway. 
* MCSS is designed to find multiple mutated driver pathways without pre-specific the number of pathways and the number of genes.
* MCSS is flexible and can borrow information from multiple types of data. We provide an example to integrate mutation data with gene expression data. The corresponding codes have been provided accordingly.


## Installation
To let researchers revise and extend MCSS more easily, we deliberately provide the source codes instead of a toolbox. Researchers can either copy the source files into their current working directory or use the following codes to add a path. The source file is fully tested in matlab 2012a and may fail in other versions of the matlab. Please send us an email (wuxx0845@umn.edu) when you meet any problems.

```matlab
addpath(genpath(<path to your MCSS directory>))
```

## MCSS

### Some background

Based on the corresponding new concepts of coverage and mutual exclusivity, MCSS is designed for de novo discovery of mutated driver pathways in cancer. Since the computational problem is a combinatorial optimization with an objective function involving a discontinuous indicator function in high dimension, many existing optimization algorithms, such as a brute force enumeration, gradient descent and Newton’s methods, are practically infeasible or directly inapplicable. That's why we develop MCSS based on a novel formulation of the problem as non-convex programming and non-convex regularization. The method is computationally more efficient, effective and scalable than existing Monte Carlo searching and several other algorithms. For more details, see our upcoming paper.

* Binghui Liu, Chong Wu, Xiaotong Shen, and Wei Pan. "CA novel and efficient algorithm for de novo discovery of mutated driver pathways in cancer." Under revision.

### How to use

* First, we generated a simple n by p mutation data set A with a 1 indicating a mutation and 0 otherwise, where n and p is the number of observations and the number of genes, respectively.  For each patient, a gene in a driver pathway was randomly selected and it mutated with probability p1, and another gene in the driver pathway was randomly selected to have a mutation probability p2. Other genes outside the pathway mutated with probability p3. We hope the following simple simulation codes could help researchers do some further research.

```matlab
n = 100; p = 500; p1 = 0.95; p2 = 0.01; p3 = 0.05

A=zeros(n,p); % A is n by p matrix

s=RandStream('mcg16807','Seed',1); % set the seed.
RandStream.setGlobalStream(s);
    
M1=1:4;
beta_M1=zeros(p,1);
beta_M1(M1,1)=ones(length(M1),1);
for i=1:n
	g1=randperm(length(M1));
	g1=M1(g1(1));
	A(i,g1)=random('bino',1,p1,1,1);
	if A(i,g1)==1
		M2=setdiff(M1,g1);
		g2=randperm(length(M1)-1);
		g2=M2(g2(1));
		A(i,g2)=random('bino',1,p2,1,1);
    end
end

M3=setdiff(1:p,M1);
for j=1:(p-length(M1))
	A(:,M3(j))=random('bino',1,p3,n,1);
end

A(1:5,1:5)
```

* Second, we apply the main function *MCSS* to get the estimates. The non-zero estimates stand for the driver mutations.

```matlab
beta_new = MCSS(A,1,1,0.1);

%	function beta_new = MCSS(A,lambda,tau1,tau2,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
% Input: 
%		A--n*p mutation matrix
%		lambda--a tuning parameter controlling the sparseness of the solutions
%		tau1--a tuning parameter in truncated lasso penalty
%		tau2--a tuning parameter relates to lambda
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
%
% Output: 
%		beta_new--a resulting estimate of beta.
%
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0
```


## Selecting the tuning parameters 
Selecting the tuning parameters is an important step for using MCSS. A better choice of tuning parameters will lead to a better result. Here, we provide a simple function that can select the tuning parameters based on the cross validation automatically. Since the results only depend on the ratio of tau1/tau2, we only search tau2 and fix tau1.


```matlab
Lambda=(1:4:10)*1;

tau1=1;

tau2=0.1;
[beta_hat,M_list,CF_list] =  MCSS_CV1(A,Lambda,tau2,tau1);
%function [beta_t,CF] = MCSS_CV1(A,Lambda,Tau2, tau1,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
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
%
% Output: 
%		beta_t--a resulting estimate of beta
%       CF--cost values corresponding Lambda
%
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0

```

## Searching solutions with multiple starting values

MCSS is able to find multiple mutated driver pathways at the same time. Unlike existing methods, which usually need pre-specific the number of pathways and the number of the genes in each pathway, MCSS can automatically find multiple mutated driver pathways which have low cost. The idea is that we can apply MCSS with different starting values and hopefully we can find multiple solutions that have low cost. We provide some codes to generate a series of random starting values automatically and then apply *MCSS_CV2* to select the "best" tuning parameters. Based on our empirical experience, using solutions from a subset of genes that have relatively large mutation rates may yield better solutions. The input of *MCSS_CV2* is almost the same as *MCSS_CV1* except beta_0 is a p by q matrix which contains q different starting points.

```matlab
Lambda=(1:4:10)*1;

tau1=1;

tau2=0.1;
alpha = 0.001;

%% Generate random starting values
beta_0_start = zeros(p,T)
for tt = 1:T
	randIndex = tt;
    s=RandStream('mcg16807','Seed',randIndex);
    RandStream.setGlobalStream(s);
    
    beta_0=random('bino',1,randi([2,8])/p,p,1);
    b5=find(beta_0>0)';
    while(ismember(sum(log([b5,0.5])),b6)==1)
        beta_0=random('bino',1,randi([2,8])/p,p,1);
        b5=find(beta_0>0)';
    end

    beta_0=zeros(p,1);beta_0(b5)=1;
    beta_0_start(:,tt) = beta_0;
end

[beta_hat,M_list,CF_list] = MCSS_CV2(A,Lambda,tau2,tau1,alpha, beta_0_final);

% function  [beta_hat,M_list,CF_list] = MCSS_CV2(A,Lambda,Tau2, tau1,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
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
%
% Output: 
%		beta_t--a resulting estimate of beta
%       CF--cost values corresponding Lambda
%
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0
```


## Mutation + gene expression data

Since MCSS is under the penalization regression framework, other types of genomic data, such as gene expression data, can be easily included. Here, we show a simple example and provide the corresponding codes（*MCSS_ME*, *MCSS\_ME\_CV1*, and  *MCSS\_ME\_CV2*). 

Note that the genes in the same pathway usually collaborate with each other, the expression of the genes in the same pathway usually have higher correlation than those from different pathways. Taking this information into account, we revise the object function to consider the correlation among the gene expression data. See more details in Integrative analysis subsection (2.6) in our upcoming paper.

Along the way, we developed the *MCSS_ME*, *MCSS_ME_CV1*, and  *MCSS_ME_CV2* for implementing MCSS algorithm and selecting the tuning parameters with both mutation and gene expression data.

Unlike the *MCSS*, in *MCSS_ME* we need choose a suitable *gamma* to balance the contributions to the new cost function from mutation data and from gene expression data. After determining *gamma*, we can use *MCSS_ME*, *MCSS_ME_CV1* and *MCSS_ME_CV2* almost the same as before.

### Generating mutation and gene expression data
```matlab
% generating mutation data
n = 100; p = 1000; p1 = 0.95; p2 = 0.01; p3 = 0.05

A=zeros(n,p); % A is n by p matrix

s=RandStream('mcg16807','Seed',1); % set the seed.
RandStream.setGlobalStream(s);
    
M1=1:4;
beta_M1=zeros(p,1);
beta_M1(M1,1)=ones(length(M1),1);
for i=1:n
	g1=randperm(length(M1));
	g1=M1(g1(1));
	A(i,g1)=random('bino',1,p1,1,1);
	if A(i,g1)==1
		M2=setdiff(M1,g1);
		g2=randperm(length(M1)-1);
		g2=M2(g2(1));
		A(i,g2)=random('bino',1,p2,1,1);
    end
end

M3=setdiff(1:p,M1);
for j=1:(p-length(M1))
	A(:,M3(j))=random('bino',1,p3,n,1);
end

A(1:5,1:5)

% Generating gene expression data, such that for the genes in the same set, their expression levels are highly correlated (0.9); for the genes in the different set, their expression levels are almost independent (0.1).

pG = p;
V=ones(pG,pG);
n_eig0=-1;
nnn=0;
while n_eig0 < 0 
	V=unifrnd(0.1,0.1,pG,pG); sum_random=4; 
	V(1:4,1:4)=unifrnd(0.9-0.0,0.9+0.0,4,4); 
	while pG-sum_random>20 
		a_random=randi([2,20]);
		b_random=unifrnd(0.9,0.9);
		V((sum_random+1):(sum_random+a_random),(sum_random+1):(sum_random+a_random))=unifrnd(b_random-0.0,b_random+0.0,a_random,a_random);
		sum_random=sum_random+a_random;
	end

	a_random=pG-sum_random;
	b_random=unifrnd(0.9,0.9);
	V((sum_random+1):(sum_random+a_random),(sum_random+1):(sum_random+a_random))=unifrnd(b_random-0,b_random+0,a_random,a_random);
	V=tril(V); V=V+eye(pG); V=V+V'; V=min(V,1);
	n_eig0=all(diag(eig(V)));
	nnn=nnn+1;
end    

V(1:10,1:10)
G=mvnrnd(zeros(pG,1),V,n); 
```

### MCSS_ME

*MCSS_ME* combines the information from both mutation and gene expression data and returns the beta value.

```matlab
beta_new = MCSS_ME(A,G,0.5,1,1,0.1);
%function beta_new =MCSS_ME(A,G,gamma,lambda,tau1,tau2,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
% Input:
%		A--n*p mutation matrix
%		G--n*p gene expression matrix
%		gamma--a tuning parameter controlling the contributions between mutation data and gene expression data
%		lambda--a tuning parameter controlling the sparseness of the solutions
%		tau1--a tuning parameter in truncated lasso penalty
%		tau2--a tuning parameter relates to lambda
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
%
% Output: beta_new--a resulting estimate of beta.
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0
```

### MCSS\_ME\_CV1

*MCSS\_ME\_CV1* works similar as *MCSS\_CV1*, selecting the tuning parameters based on cross validation.

```matlab
Lambda=(1:4:10)*1;

tau1=1;

tau2=0.1;
[beta_hat,M_list,CF_list] =  MCSS_ME_CV1(A,G,0.5,Lambda,tau2,tau1);

%function [beta_t,CF] = MCSS_ME_CV1(A,G, gamma, Lambda,Tau2, tau1,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
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
%
% Output: 
%		beta_t--a resulting estimate of beta
%       CF--cost values corresponding Lambda
%
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0

```

### MCSS\_ME\_CV2

*MCSS\_ME\_CV2* works similar as *MCSS\_CV2*. Based on multiple starting points, *MCSS\_ME\_CV2* finds some gene pathways with low cost.

```matlab
Lambda=(1:4:10)*1;

tau1=1;

tau2=0.1;
alpha = 0.001;

%% Generate random starting values
beta_0_start = zeros(p,T)
for tt = 1:T
	randIndex = tt;
    s=RandStream('mcg16807','Seed',randIndex);
    RandStream.setGlobalStream(s);
    
    beta_0=random('bino',1,randi([2,8])/p,p,1);
    b5=find(beta_0>0)';
    while(ismember(sum(log([b5,0.5])),b6)==1)
        beta_0=random('bino',1,randi([2,8])/p,p,1);
        b5=find(beta_0>0)';
    end

    beta_0=zeros(p,1);beta_0(b5)=1;
    beta_0_start(:,tt) = beta_0;
end

[beta_hat,M_list,CF_list] = MCSS_CV2(A,G,0.5, Lambda,tau2,tau1,alpha, beta_0_final);

% function  [beta_hat,M_list,CF_list] = MCSS_ME_CV2(A,G,gamma,Lambda,Tau2, tau1,alpha, beta_0, method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
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
%
% Output: 
%		beta_t--a resulting estimate of beta
%		M_list--resulting list of gene sets corresponding to Lambda
%		CF--cost values corresponding to Lambda
%
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0
```

## Stay In Touch

- For latest releases and announcements, follow on GitHub: [@ChongWu-Biostat](https://github.com/ChongWu-Biostat)


## License

[MIT](http://opensource.org/licenses/MIT)

Copyright (c) 2013-present, Chong Wu

