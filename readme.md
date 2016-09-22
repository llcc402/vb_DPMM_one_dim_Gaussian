# Bayesian Variational Inference for Dirichlet Process Mixture of Gaussians

## Model description

We generate one dimensional observations from a finite mixture of Gaussians. Then we approximate the posterior by a fully factorized distribution who is a product of a series of unary distributions. Our task is to compute the parameters of these distributions.

The model is build like:

        X_i|Z_i,mu,lambda    ~   N(mu_{Z_i}, lambda_{Z_i}^{-1})   i = 1,...,n
        Pr(Z_i = j|pi)       =   pi_j \prod_{t=1}^{j-1} (1 - pi_t)  j = 1,2,...
        pi_j|alpha           ~   Beta(1, alpha)   j = 1,2,...
        (mu_j, lambda_j|a,b) ~   N(mu_j|0, lambda_j) Gamma(lambda_j|a, b)  j = 1,2,...

The posterior is assumed to be:

       mu_j      ~   N(mu_j|mu(j).mean, mu(j).precision^{-1})          
       lambda_j  ~   Gamma(lambda_j|lambda(j).a, lambda(j).b) 
       pi_j      ~   Beta(pi_j|pi(j).a1, pi(j).a2)

## Functions included

1. data_generate: It has 2 parameters, the first parameter specify the number of clusters and the second parameter specify the number of observations in each cluster. We fix the mean of the whole dataset is 0. The distance between the centers of two clusters is at least 5. The standard deviation of each cluster is distributed as a Gamma distribution with mean 0.5 and variance 0.25.

2. vb: It has 6 parameters. We use the value of the mean of the first clster in each iteration to determine whether the algorithm converges. 


## Scripts included

1. main: Run this script for example

## Pictures included

1. mean_iter: The mean of the first cluster v.s. iterations.

2. data_hist: The histogram of the dataset.

3. inferred_z: The inferred indicator for each observation. The height of the points is just used to distiguish the clusters (See main.m for details). 
