# MCSS  

Cancer is characterized by numerous somatic mutations, however, only a subset of mutations, called *driver* mutations, contributes to tumor growth and progression. Minimum cost subset selection (**MCSS**) is a new method for de novo discovery of mutated driver pathway in cancer. 

## Features

* We provide a novel approximation to a combinatorial problem through regularization. We provide two versions of the difference of convex algorithms to solve the non-convex objective function.
* MCSS is designed to find multiple mutated driver pathways.
* MCSS is flexible and can take multiple type of information into account. We provide an example to integrate mutation data with gene expression data. Some other possible information, such as gene-gene interaction network may be incorporated with minor revisions. 

## Installation
To let researchers revise and extend MCSS and learn the details of MCSS, we deliberately provide the source codes instead of a toolbox. You can either copy the source files into your current working directory or use the following codes to add a path. The source file is fully tested in matlab 2012a and may fail in other versions of the matlab. Please send us an email (wuxx0845@umn.edu) when you meet any problems.

```matlab
addpath(genpath(<path to your MCSS directory>))
```

## Using MCSS

### Some background

Based on the corresponding new concepts of coverage and mutual exclusivity, MCSS designed for de novo discovery of mutated driver pathways in cancer. Since the computational problem is a combinatorial optimization with an objective function involving a discontinuous indicator function in high dimension, many existing optimization algorithms, such as a brute force enumeration, gradient descent and Newton’s methods, are practically infeasible or directly inapplicable. That's why we develop MCSS based on a novel formulation of the problem as non-convex programming and non-convex regularization. The method is computationally more efficient, effective and scalable than existing Monte Carlo searching and several other algorithms. For more details, see our upcoming paper.

* Binghui Liu, Chong Wu, Xiaotong Shen, and Wei Pan. "CA novel and e cient algorithm for de novo discovery of mutated driver pathways in cancer." In revision.

### How to use

* First, we want to generate a simple n by p mutation data set A with a 1 indicating a mutation and 0 otherwise, where n is the number of observations, p is the number of genes. We provide the following simple simulation codes that hopefully can help you do some further research. For each patient, a gene in a driver pathway was randomly selected and it mutated with probability p1, and another gene in the driver pathway was randomly selected to have a mutation probability p2. Other genes outside the pathway mutated with probability p3.


```matlab

```

* Second, we can use the main function *MCSS* to get the estimates. The non-zero estimates correspond to the driver mutations.

```matlab
beta_new = MCSS(A,1,1,0.1);

%beta_new = MCSS(A,lambda,tau1,tau2,alpha, beta_0, 
%			method,weight,max_iter_num, Err,max_iter_num_s,Err_s)
% Input: 
%		A--n*p mutation matrix
%		lambda--tuning parameters controlling the sparseness of the solutions
%		tau1--tuning parameters in truncated lasso penalty
%		tau2--tuning parameters relates to lambda
%		alpha--tuning parameter with ridge penalty
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
```


## Selecting the tuning parameters 
Selecting the tuning parameters is an important step for using MCSS. A better choice of tuning parameters will lead to a better result. Here, we provide a simple function that can help you select the tuning parameters based on the cross validation automatically. Note that the results only depends on the ratio of tau1/tau2, thus why we only select the tau2 and put tau1 as a fixed value.


```matlab
Lambda=(1:4:10)*1;

tau1=1;

tau2=0.1;
MCSS_CV1(A,Lambda,tau2,tau1)
%function [beta_t,CF] = MCSS_CV1(A,Lambda, Tau2 ,tau1
%							,alpha, beta_0, method,weight,max_iter_num
%							, Err,max_iter_num_s,Err_s)
% Input:
%		A--n*p mutation matrix
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
% Output: beta_t--a resulting estimate of beta
%         CF--cost values corresponding Lambda
% Author: Binghui Liu and Chong Wu (wuxx0845@umn.edu)
% Maintainer: Chong Wu (wuxx0845@umn.edu)
% Version: 1.0
```

## Search solutions with multiple starting values

One feature of the MCSS is that it is able to find multiple mutated driver pathways at the same time. Unlike existing methods, which usually need pre-specific the number of pathways and the number of the genes in each pathway beforehand, MCSS can automatically find multiple mutated driver pathways which have low cost. The idea is that we can apply MCSS with different starting values and hopefully we can find multiple solutions that have low cost. We provide a function to generate a series of random starting values automatically and then apply the cross validation to select the "best" tuning parameters. The input of *MCSS_CV2* is exactly the same as *MCSS_CV1*.


```matlab
Lambda=(1:4:10)*1;

tau1=1;

tau2=0.1;
MCSS_CV2(A,Lambda,tau2,tau1)
```


## Integrate other data sets
Since MCSS is under the penalization regression framework, other type of genomic data, such as gene expression data, can be easily included. Here, we just show a simple example and provide the corresponding codes （*MCSS_ME*). You can learn more details from our upcoming paper. 



